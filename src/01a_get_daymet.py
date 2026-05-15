import pydaymet as daymet
import pandas as pd
from tqdm import tqdm
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
import time

# -------------------------
# Configuration
# -------------------------
MAX_WORKERS = 4       # Reduce if API errors increase
SLEEP_TIME = 0     # Seconds between requests, increase if rate limited
INPUT_PATH = "data/raw/covariates/daymet/daymet_sites.csv"
OUTPUT_PATH = "data/raw/covariates/daymet/daymet_output.csv"
FAILED_PATH = "data/raw/covariates/daymet/daymet_failed.csv"
SAVE_INTERVAL = 100   # Save progress every N successful results
VARIABLES = ["tmin", "vp", "prcp", "dayl"]
# -------------------------

# Read in sites and parse date column
sites = pd.read_csv(INPUT_PATH)
sites["night"] = pd.to_datetime(sites["night"])

# Deduplicate to avoid redundant API calls
unique_sites = sites.drop_duplicates(subset=["latitude", "longitude", "night"])

results = []
failed = []

def fetch_daymet(row):
    time.sleep(SLEEP_TIME)
    data = daymet.get_bycoords(
        (row["longitude"], row["latitude"]),
        dates=(row["night"], row["night"]),
        variables=VARIABLES
    )
    # Filter to just the deployment date
    data = data[data.index == row["night"]]
    data["location_name"] = row["location_name"]
    data["night"] = row["night"]
    return data

with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
    futures = {executor.submit(fetch_daymet, row): row for _, row in unique_sites.iterrows()}
    
    for future in tqdm(as_completed(futures), total=len(futures)):
        row = futures[future]
        try:
            data = future.result()
            results.append(data)

            if len(results) % SAVE_INTERVAL == 0:
                pd.concat(results).to_csv(OUTPUT_PATH, index=False)

        except Exception as e:
            failed.append({
                "location_name": row["location_name"],
                "night": row["night"],
                "latitude": row["latitude"],
                "longitude": row["longitude"],
                "error": str(e),
                "timestamp": datetime.now()
            })
            pd.DataFrame(failed).to_csv(FAILED_PATH, index=False)

# Final save
daymet_out = pd.concat(results)
daymet_out.to_csv(OUTPUT_PATH, index=False)
print(f"Done! Retrieved data for {len(daymet_out)} sites.")
print(f"Failed for {len(failed)} sites. See {FAILED_PATH} for details.")