import rasterio
from rasterio.warp import transform, CRS

tif_path = r"C:\Users\Mias - Laptop\Downloads\XFOIL6.99\NAM_combined-Weibull-k_50m.tif"

lat, lon = -26.6453, 14.1604

with rasterio.open(tif_path) as ds:

    src_crs = ds.crs or CRS.from_epsg(4326)  # <-- FIX HERE

    x, y = transform(
        "EPSG:4326",
        src_crs,
        [lon],
        [lat]
    )

    k = list(ds.sample([(x[0], y[0])]))[0][0]

print("Weibull k:", k)