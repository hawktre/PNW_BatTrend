#!/usr/bin/env python3
import os
import time
import pandas as pd
import pydaymet as dm

INPUT_PATH = "data/raw/covariates/daymet/daymet_sites.csv"
OUTPUT_PATH = "data/raw/covariates/daymet/daymet_output.csv"
CHUNK_SIZE = 200

def process_chunk_data(clm, loc_names_chunk, coords_chunk, year_df, new_records):
    # Normalize index to timezone-naive date
    clm.index = pd.to_datetime(clm.index).normalize()
    
    for idx, loc_name in enumerate(loc_names_chunk):
        # Find all nights we need for this location in the current year
        loc_nights = year_df[year_df["location_name"] == loc_name]
        
        for _, row in loc_nights.iterrows():
            night_dt = pd.to_datetime(row["night"]).normalize()
            if night_dt in clm.index:
                tmin = clm.loc[night_dt, (idx, "tmin (degrees C)")]
                vp = clm.loc[night_dt, (idx, "vp (Pa)")]
                prcp = clm.loc[night_dt, (idx, "prcp (mm/day)")]
                dayl = clm.loc[night_dt, (idx, "dayl (s)")]
                
                new_records.append({
                    "location_name": loc_name,
                    "night": night_dt.strftime("%Y-%m-%d"),
                    "latitude": row["latitude"],
                    "longitude": row["longitude"],
                    "tmin": float(tmin),
                    "vp": float(vp),
                    "prcp": float(prcp),
                    "dayl": float(dayl)
                })

def main():
    print("Starting Daymet data retrieval using pydaymet...")
    
    # 1. Load input sites
    if not os.path.exists(INPUT_PATH):
        raise FileNotFoundError(f"Input file not found at {INPUT_PATH}")
        
    sites_df = pd.read_csv(INPUT_PATH)
    sites_df = sites_df.dropna(subset=["latitude", "longitude", "night"])
    sites_df["night"] = pd.to_datetime(sites_df["night"])
    sites_df["year"] = sites_df["night"].dt.year
    
    # 2. Check for completed records to skip (differential downloading)
    completed_keys = set()
    completed_df_list = []
    
    if os.path.exists(OUTPUT_PATH):
        try:
            existing_output = pd.read_csv(OUTPUT_PATH)
            # Find rows that are fully complete (all weather variables are non-null)
            required_cols = ["tmin", "vp", "prcp", "dayl"]
            if all(col in existing_output.columns for col in ["location_name", "night"] + required_cols):
                complete_rows = existing_output.dropna(subset=required_cols)
                for _, row in complete_rows.iterrows():
                    # Format date to match
                    d_str = pd.to_datetime(row["night"]).strftime("%Y-%m-%d")
                    completed_keys.add((row["location_name"], d_str))
                    
                completed_df_list.append(complete_rows)
                print(f"Found {len(complete_rows)} completed site-nights in existing output. Skipping them.")
        except Exception as e:
            print(f"Could not read existing output file: {e}. Starting fresh.")
            
    # Filter for pending records
    pending_rows = []
    for _, row in sites_df.iterrows():
        d_str = row["night"].strftime("%Y-%m-%d")
        if (row["location_name"], d_str) not in completed_keys:
            pending_rows.append(row)
            
    if not pending_rows:
        print("All records already downloaded. Nothing to do!")
        return
        
    pending_df = pd.DataFrame(pending_rows)
    print(f"Found {len(pending_df)} pending site-night combinations to download.")
    
    # 3. Group pending records by year and retrieve in chunks
    new_records = []
    years = sorted(pending_df["year"].unique())
    
    for year in years:
        year_df = pending_df[pending_df["year"] == year]
        
        # Get unique coordinates for this year to minimize API calls
        unique_coords_df = year_df.drop_duplicates(subset=["location_name", "latitude", "longitude"])
        coords_all = list(zip(unique_coords_df["longitude"], unique_coords_df["latitude"]))
        loc_names_all = unique_coords_df["location_name"].tolist()
        
        n_unique = len(coords_all)
        print(f"Year {year}: Querying {n_unique} unique coordinates...")
        
        # Query in chunks of CHUNK_SIZE to be safe and avoid HTTP timeouts
        for start_idx in range(0, n_unique, CHUNK_SIZE):
            end_idx = min(start_idx + CHUNK_SIZE, n_unique)
            coords_chunk = coords_all[start_idx:end_idx]
            loc_names_chunk = loc_names_all[start_idx:end_idx]
            
            n_chunk = len(coords_chunk)
            print(f"  Retrieving chunk [{start_idx + 1}-{end_idx}] (size: {n_chunk})...")
            
            try:
                start_time = time.time()
                clm = dm.get_bycoords(
                    coords=coords_chunk,
                    dates=int(year),
                    variables=["tmin", "vp", "prcp", "dayl"],
                    validate_filesize=False  # Speeds up cached lookups
                )
                print(f"  Chunk download completed in {time.time() - start_time:.2f} seconds.")
                
                # If only 1 coordinate, pydaymet returns single Index columns, so we promote to MultiIndex
                if n_chunk == 1:
                    clm.columns = pd.MultiIndex.from_product([[0], clm.columns])
                    
                process_chunk_data(clm, loc_names_chunk, coords_chunk, year_df, new_records)
                
            except Exception as e:
                print(f"  Error downloading chunk: {e}. Retrying coordinates individually...")
                for idx, (lon, lat) in enumerate(coords_chunk):
                    loc_name = loc_names_chunk[idx]
                    try:
                        single_clm = dm.get_bycoords(
                            coords=(lon, lat),
                            dates=int(year),
                            variables=["tmin", "vp", "prcp", "dayl"],
                            validate_filesize=False
                        )
                        single_clm.columns = pd.MultiIndex.from_product([[0], single_clm.columns])
                        process_chunk_data(single_clm, [loc_name], [(lon, lat)], year_df, new_records)
                    except Exception as e_single:
                        print(f"    Failed for site {loc_name} ({lat}, {lon}): {e_single}")

    # 4. Save and compile final results
    if new_records:
        new_df = pd.DataFrame(new_records)
        completed_df_list.append(new_df)
        
    if completed_df_list:
        final_df = pd.concat(completed_df_list, ignore_index=True)
        # Deduplicate just in case
        final_df = final_df.drop_duplicates(subset=["location_name", "night"])
        # Sort by location_name and night for clean representation
        final_df = final_df.sort_values(by=["location_name", "night"])
        
        # Ensure parent directories exist
        os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
        final_df.to_csv(OUTPUT_PATH, index=False)
        print(f"Done! Retrieved data saved to {OUTPUT_PATH}")
        print(f"Total compiled records: {len(final_df)}")
    else:
        print("No records compiled.")

if __name__ == "__main__":
    main()
