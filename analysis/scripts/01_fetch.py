"""Download the ERA5 subset from ARCO and assemble data/era5_subset.nc."""

from arbudget import fetch

if __name__ == "__main__":
    fetch.fetch(workers=6)
    ds = fetch.assemble()
    print(ds)
