# %% load MATLAB -> cache -> summary  (run ONCE per data update)
import mi_io, config, importlib, pandas as pd
importlib.reload(config)
# %% 
NEURONS, REACH_PX = mi_io.parse_mat(config.DATA_MAT, config.struct_names())
cache = config.cache_path()
mi_io.save(NEURONS, REACH_PX, cache)
print(f"\nloaded  {config.DATA_MAT}")
print(f"cache   {cache}")
print(f"neurons {len(NEURONS)}  |  reach_px sessions {len(REACH_PX)}\n")
# %% 
s = pd.DataFrame([{k: n[k] for k in ("region", "mouse", "sess", "n_trials")} for n in NEURONS])
print(s.groupby("region").agg(n_neurons=("n_trials","size"),
                              sessions=("sess","nunique"),
                              trials_median=("n_trials","median")))
print("\nper session:")
print(s.groupby(["region","mouse","sess"]).size().rename("n_neu").reset_index())

