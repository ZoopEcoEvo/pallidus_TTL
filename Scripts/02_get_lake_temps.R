# Download continuous water temperature for Cheney Reservoir, Kansas
# (USGS 07144790, "Cheney Reservoir near Cheney, KS") from the USGS Water Data
# APIs and save it to Raw_data/lake_temps/.
#
# Run this once (or when you want newer data) from the project root, or via
# 00_controller.R with download_lake_temps = T. The report reads the saved CSV,
# so it does not need internet access.
#
# Requirements
#   - dataRetrieval >= 2.7.22 (read_waterdata_continuous() was added in 2.7.22).
#     The package is pure R, so install.packages("dataRetrieval") works even on
#     older R versions (it may install from source).
#   - Recommended: a free USGS API token, to avoid "HTTP 429 Too Many Requests".
#     Request one (free, emailed instantly) at https://api.waterdata.usgs.gov/signup/, then add
#     API_USGS_PAT=your-token to ~/.Renviron (usethis::edit_r_environ()) and
#     restart R.
#
# Notes
#   - Uses the modern Water Data APIs. The legacy readNWISuv() / WaterServices
#     endpoints are being retired in early 2027.
#   - The service returns at most 3 years of continuous data per request, so the
#     download runs one calendar year at a time.
#   - Times are returned in UTC and saved in UTC; convert to local time
#     (America/Chicago) when analysing.

library(dataRetrieval)
library(tidyverse)

stopifnot("dataRetrieval >= 2.7.22 is required (read_waterdata_continuous)" =
            packageVersion("dataRetrieval") >= "2.7.22")

# Settings -------------------------------------------------------------------
site      = "USGS-07144790"   # Cheney Reservoir near Cheney, KS
pcode     = "00010"           # water temperature, deg C
years     = NULL              # NULL = full period of record; or e.g. 2015:2025
ts_id     = "8dbbfd775ce04a8aaae1f9118002c32d" # set to a time_series_id if the site has more than one temperature series
out_dir   = "Raw_data/lake_temps"
out_file  = file.path(out_dir, "cheney_07144790_water_temp.csv")
meta_file = file.path(out_dir, "cheney_07144790_ts_metadata.csv")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Which continuous temperature series exist at the site? -------------------
# A site can have more than one water temperature series (e.g. sensors at
# different depths). The metadata usually describes each one; check it before
# choosing.
meta = read_waterdata_ts_meta(monitoring_location_id = site,
                              parameter_code = pcode)
if (inherits(meta, "sf")) meta = sf::st_drop_geometry(meta)
meta = as_tibble(meta)

# dataRetrieval renames the service's `id` field to `time_series_id`; accept either
id_col = intersect(c("time_series_id", "id"), names(meta))[1]
if (is.na(id_col)) stop("Could not find a time series id column in the metadata.")
meta = rename(meta, series_id = all_of(id_col))

# keep the instantaneous ("Points") series, i.e. not daily statistics
if ("computation_period_identifier" %in% names(meta)) {
  meta = filter(meta, computation_period_identifier == "Points")
}
if (nrow(meta) == 0) stop("No continuous water temperature series found for ", site)

meta_print = meta |> select(any_of(c("series_id", "parameter_name", "begin", "end",
                                     "sublocation_identifier", "web_description",
                                     "statistic_id", "unit_of_measure")))
message("Continuous water temperature series at ", site, ":")
print(meta_print, width = Inf)
write_csv(meta |> mutate(across(where(is.list), \(x) map_chr(x, \(v) paste(v, collapse = "; ")))),
          meta_file)

if (is.null(ts_id)) {
  if (nrow(meta) > 1) {
    stop("This site has ", nrow(meta), " continuous temperature series (see above and ", meta_file, "). ",
         "Set `ts_id` to the one you want (e.g. the sensor depth closest to where S. pallidus lives) and rerun.")
  }
  ts_id = meta$series_id[1]
}
series = filter(meta, series_id == ts_id)
if (nrow(series) != 1) stop("ts_id '", ts_id, "' not found in the site metadata.")

# 2. Download one calendar year at a time --------------------------------------
if (is.null(years)) {
  years = seq(year(as.Date(series$begin)), year(if (is.na(series$end)) Sys.Date() else as.Date(series$end)))
}

get_year = function(y) {
  message("Downloading ", y, " ...")
  out = tryCatch(
    read_waterdata_continuous(
      monitoring_location_id = site,
      parameter_code = pcode,
      time_series_id = ts_id,
      time = c(paste0(y, "-01-01T00:00:00Z"), paste0(y, "-12-31T23:59:59Z"))),
    error = function(e) {
      warning("Year ", y, " failed: ", conditionMessage(e), call. = FALSE)
      NULL
    })
  Sys.sleep(1)   # be gentle with the API
  out
}

raw = map(years, get_year) |> compact() |> list_rbind()
if (is.null(raw) || nrow(raw) == 0) stop("No data returned. Check the site, series and your API token.")

# 3. Tidy and save ------------------------------------------------------------
# qualifier can come back as a list column; flatten it to text so it saves cleanly
flatten = \(x) if (is.list(x)) map_chr(x, \(v) paste(na.omit(v), collapse = "; ")) else as.character(x)

lake_temps = raw |>
  as_tibble() |>
  select(any_of(c("time_series_id", "monitoring_location_id", "parameter_code",
                  "time", "value", "unit_of_measure", "approval_status", "qualifier"))) |>
  mutate(qualifier = flatten(qualifier),
         time = with_tz(as_datetime(time), "UTC")) |>
  distinct(time, .keep_all = TRUE) |>
  arrange(time)

write_csv(lake_temps, out_file)

message("Saved ", nrow(lake_temps), " readings (", format(min(lake_temps$time)), " to ",
        format(max(lake_temps$time)), " UTC) to ", out_file)
message("Sampling interval (minutes): ",
        paste(names(sort(table(as.numeric(diff(lake_temps$time), units = "mins")), decreasing = TRUE))[1:3],
              collapse = ", "), " (most common first)")
