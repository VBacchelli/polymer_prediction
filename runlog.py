"""Salva su disco configurazione e risultati delle esecuzioni.

Ogni esecuzione crea una cartella dentro `runs/`:

    runs/2026-07-14_1503__biceranoPolymers__md/
        run.json           configurazione, dati e risultati completi
        report.md          riepilogo leggibile del run
        tables/*.csv       tabelle prodotte durante l'analisi
        figures/*.png      figure salvate dal notebook

Inoltre aggiunge una riga a `runs/index.csv`, cosi' i run si possono confrontare
senza dover aprire ogni singolo file JSON.

Il nome della cartella riassume la configurazione usata:
`<timestamp>__<collection mongo>__<md|nomd>__<tag>`.

Uso tipico:

    log = RunLogger(config={...}, tag="prova")   # crea la cartella del run
    log.dataset(n_polymers=139, ...)             # registra la sorgente dei dati
    log.result("nested_cv", {"MAE_mean": 24.8})  # salva una sezione di risultati
    log.table("nested_cv_folds", nested_df)      # salva una tabella
    log.figure("shap_summary")                   # salva la figura corrente
    log.finish()                                 # chiude il run e aggiorna l'indice

Il file `run.json` viene aggiornato a ogni chiamata. Se il notebook si interrompe,
i risultati gia' prodotti restano disponibili.

`MONGO_URI` non viene mai scritto nei log perche' contiene le credenziali.
Vengono salvati solo il nome del database e quello della collection.
"""

import json
import platform
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_NAME = "polymer_prediction"


def _find_project_root(start):
    """Risale le directory fino a trovare la root del progetto."""
    start = Path(start).resolve()
    for directory in (start, *start.parents):
        if directory.name == PROJECT_NAME:
            return directory
    raise RuntimeError(
        f"Directory root {PROJECT_NAME!r} non trovata partendo da {start}"
    )


REPO_ROOT = _find_project_root(Path(__file__).parent)
RUNS_ROOT = REPO_ROOT / "runs"


def _jsonable(o):
    """Converte in valori JSON i tipi piu' comuni di NumPy e pandas."""
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, (np.bool_, bool)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, pd.Series):
        return o.to_dict()
    if isinstance(o, pd.DataFrame):
        return o.to_dict(orient="records")
    if isinstance(o, (datetime, pd.Timestamp)):
        return o.isoformat()
    if isinstance(o, Path):
        return str(o)
    return str(o)


def _git_commit():
    """Restituisce il commit usato per produrre il run, se disponibile."""
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=REPO_ROOT,
            capture_output=True, text=True, timeout=5, check=True,
        )
        return out.stdout.strip()
    except Exception:
        return None


def _pkg_versions():
    versions = {"python": platform.python_version()}
    for name in ["numpy", "pandas", "scikit-learn", "rdkit", "shap", "umap-learn"]:
        try:
            from importlib.metadata import version

            versions[name] = version(name)
        except Exception:
            pass
    return versions


def _slug(text):
    keep = [c if (c.isalnum() or c in "-_") else "-" for c in str(text)]
    return "".join(keep).strip("-")


class RunLogger:
    """Raccoglie e salva i dati principali di un'esecuzione."""

    def __init__(self, config, tag=None, root=RUNS_ROOT):
        cfg = dict(config)
        stamp = datetime.now().strftime("%m-%d_%H%M")
        # Il token indica il tipo di dati MD usati; "nomd" indica che non sono state usate.
        md_token = _slug(cfg.get("md_dataset") or "md") if cfg.get("use_md") else "nomd"
        parts = [
            stamp,
            _slug(cfg.get("collection", "nodb")),
            md_token,
        ]
        if cfg.get("use_mordred"):
            parts.append("mrd")
        if cfg.get("md_only"):
            parts.append("mdonly")
        if tag:
            parts.append(_slug(tag))
        self.run_id = "_".join(parts)

        self.dir = Path(root) / self.run_id
        (self.dir / "tables").mkdir(parents=True, exist_ok=True)
        (self.dir / "figures").mkdir(parents=True, exist_ok=True)

        self._t0 = time.time()
        self.data = {
            "run_id": self.run_id,
            "tag": tag,
            "started_at": datetime.now().isoformat(timespec="seconds"),
            "finished_at": None,
            "duration_min": None,
            "git_commit": _git_commit(),
            "versions": _pkg_versions(),
            "config": cfg,
            "dataset": {},
            "results": {},
            "tables": {},
            "figures": {},
            "notes": [],
        }
        self._save()
        print(f"[runlog] {self.dir}")

    # -- salvataggio -------------------------------------------------------

    def _save(self):
        with open(self.dir / "run.json", "w", encoding="utf-8") as f:
            json.dump(self.data, f, indent=2, default=_jsonable, ensure_ascii=False)

    def config(self, **fields):
        """Aggiunge altri parametri alla configurazione del run."""
        self.data["config"].update(fields)
        self._save()
        return self

    def dataset(self, **fields):
        """Registra sorgente, dimensione e intervallo dei dati usati."""
        self.data["dataset"].update(fields)
        self._save()
        return self

    def result(self, section, payload):
        """Aggiunge risultati a una sezione del manifest."""
        self.data["results"].setdefault(section, {}).update(payload)
        self._save()
        return payload

    def table(self, name, df, index=True):
        """Salva una tabella CSV e ne registra il percorso nel manifest."""
        path = self.dir / "tables" / f"{name}.csv"
        df.to_csv(path, index=index)
        self.data["tables"][name] = str(path.relative_to(self.dir)).replace("\\", "/")
        self._save()
        return df

    def figure(self, name, fig=None, dpi=150):
        """Salva una figura PNG e ne registra il percorso nel manifest."""
        fig = fig or plt.gcf()
        path = self.dir / "figures" / f"{name}.png"
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        self.data["figures"][name] = str(path.relative_to(self.dir)).replace("\\", "/")
        self._save()
        return fig

    def note(self, text):
        """Aggiunge una nota libera al report finale."""
        self.data["notes"].append(str(text))
        self._save()
        return self

    # -- chiusura ----------------------------------------------------------

    def finish(self, quiet=False):
        """Chiude il run, scrive il report e aggiorna l'indice."""
        self.data["finished_at"] = datetime.now().isoformat(timespec="seconds")
        self.data["duration_min"] = round((time.time() - self._t0) / 60, 2)
        self._save()

        report = self._report()
        with open(self.dir / "report.md", "w", encoding="utf-8") as f:
            f.write(report)
        self._append_index()

        if not quiet:
            print(report)
            print(f"\n[runlog] scritto in {self.dir}")
        return self.dir

    # -- report ------------------------------------------------------------

    @staticmethod
    def _render(value, indent=0):
        """Converte un valore annidato in un elenco Markdown."""
        pad = "  " * indent
        lines = []
        if isinstance(value, dict):
            for k, v in value.items():
                if isinstance(v, (dict, list)):
                    lines.append(f"{pad}- **{k}**:")
                    lines += RunLogger._render(v, indent + 1)
                else:
                    lines.append(f"{pad}- **{k}**: {RunLogger._fmt(v)}")
        elif isinstance(value, list):
            if all(not isinstance(v, (dict, list)) for v in value):
                lines.append(f"{pad}- {', '.join(RunLogger._fmt(v) for v in value) or '(vuoto)'}")
            else:
                for v in value:
                    lines += RunLogger._render(v, indent)
        else:
            lines.append(f"{pad}- {RunLogger._fmt(value)}")
        return lines

    @staticmethod
    def _fmt(v):
        if isinstance(v, float):
            return f"{v:.4g}"
        return str(v)

    def _report(self):
        d = self.data
        ds = d["dataset"]
        cfg = d["config"]
        md_on = bool(cfg.get("use_md"))

        L = [
            f"# Run `{d['run_id']}`",
            "",
            f"- **inizio**: {d['started_at']}  |  **fine**: {d['finished_at']}"
            f"  |  **durata**: {d['duration_min']} min",
            f"- **commit**: `{d['git_commit']}`",
        ]
        if d["tag"]:
            L.append(f"- **tag**: {d['tag']}")

        L += [
            "",
            "## Dati",
            "",
            f"- **sorgente**: MongoDB `{ds.get('db', '?')}.{ds.get('collection', '?')}`",
            f"- **target**: {cfg.get('target')}",
            f"- **polimeri usati**: {ds.get('n_polymers', '?')}"
            + (f" (su {ds['n_docs']} documenti nella collection)" if ds.get("n_docs") else ""),
            "- **descrittori del monomero**: RDKit + Morgan + MACCS"
            + (" + Mordred" if cfg.get("use_mordred") else " (`USE_MORDRED = False`)"),
        ]
        if ds.get("target_min") is not None:
            L.append(f"- **range {cfg.get('target')}**: {ds['target_min']:.0f} - {ds['target_max']:.0f} K"
                     f" (media {ds['target_mean']:.0f}, std {ds['target_std']:.0f})")

        L += ["", "## Feature MD", ""]
        if md_on:
            L += [
                f"- **usate**: SI ({len(cfg.get('md_feats', []))} feature)",
                f"- **sorgente**: `{ds.get('md_source', cfg.get('md_source', '?'))}`",
                "- **lista**:",
            ]
            L += [f"  - `{f}`" for f in cfg.get("md_feats", [])]
        else:
            L += [
                "- **usate**: NO (`USE_MD = False`)",
                "- il run e' la sola pipeline RDKit: un solo braccio (`base`), niente confronto base/full.",
            ]

        L += ["", "## Configurazione", ""]
        L += self._render({k: v for k, v in cfg.items()
                           if k not in ("md_feats", "collection", "db", "target", "use_md")})

        L += ["", "## Risultati", ""]
        if not d["results"]:
            L.append("_(nessun risultato registrato)_")
        for section, payload in d["results"].items():
            L += [f"### {section}", ""]
            L += self._render(payload)
            L.append("")

        if d["tables"]:
            L += ["## Tabelle", ""]
            L += [f"- [`{v}`]({v})" for v in d["tables"].values()]
            L.append("")
        if d["figures"]:
            L += ["## Figure", ""]
            L += [f"- [`{v}`]({v})" for v in d["figures"].values()]
            L.append("")
        if d["notes"]:
            L += ["## Note", ""]
            L += [f"- {n}" for n in d["notes"]]
            L.append("")

        L += ["## Ambiente", ""]
        L += self._render(d["versions"])
        return "\n".join(L) + "\n"

    # -- indice dei run ----------------------------------------------------

    def _append_index(self):
        """Aggiunge all'indice una riga con le informazioni principali del run."""
        d, cfg, res = self.data, self.data["config"], self.data["results"]

        def g(section, key, default=None):
            return res.get(section, {}).get(key, default)

        row = {
            "run_id": d["run_id"],
            "started_at": d["started_at"],
            "duration_min": d["duration_min"],
            "tag": d["tag"],
            "git_commit": d["git_commit"],
            "db": cfg.get("db"),
            "collection": cfg.get("collection"),
            "target": cfg.get("target"),
            "n_polymers": d["dataset"].get("n_polymers"),
            "use_md": bool(cfg.get("use_md")),
            "md_only": bool(cfg.get("md_only")),
            "n_md_feats": len(cfg.get("md_feats", [])),
            "md_feats": " ".join(cfg.get("md_feats", [])),
            "md_source": d["dataset"].get("md_source"),
            "use_mordred": bool(cfg.get("use_mordred")),
            "n_mordred_feats": g("mordred", "n_tenuti"),
            "random_state": cfg.get("random_state"),
            # La metrica principale viene dalla nested CV.
            "nested_MAE_base": g("nested_cv", "md_only_MAE_mean" if cfg.get("md_only") else "base_MAE_mean"),
            "nested_MAE_base_std": g("nested_cv", "md_only_MAE_std" if cfg.get("md_only") else "base_MAE_std"),
            "nested_R2_base": g("nested_cv", "md_only_R2_mean" if cfg.get("md_only") else "base_R2_mean"),
            "nested_MAE_full": g("nested_cv", "full_MAE_mean"),
            "nested_MAE_full_std": g("nested_cv", "full_MAE_std"),
            "nested_R2_full": g("nested_cv", "full_R2_mean"),
            "delta_MAE_full_base": g("nested_cv", "delta_MAE_mean"),
            "delta_MAE_pvalue": g("nested_cv", "paired_ttest_p"),
            # Metriche del singolo split e informazioni riassuntive.
            "split_test_MAE_base": g("single_split", "md_only_test_mae" if cfg.get("md_only") else "base_test_mae"),
            "split_test_MAE_full": g("single_split", "full_test_mae"),
            "n_final_feats": g("final_model", "n_features"),
            "selection_jaccard": g("selection_stability", "jaccard_mean"),
            "md_shap_share_pct": g("md_importance", "shap_share_md_pct"),
        }
        
        index_path = Path(self.dir).parent / "index.csv"
        df = pd.DataFrame([row])
        if index_path.exists():
            # Se il formato dell'indice cambia, ricarica le righe esistenti e
            # riscrive il file usando l'unione delle colonne.
            df = pd.concat([pd.read_csv(index_path), df], ignore_index=True)
        df.to_csv(index_path, index=False)
        return index_path


def load_runs(root=RUNS_ROOT):
    """Carica l'indice dei run in un DataFrame."""
    path = Path(root) / "index.csv"
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)
