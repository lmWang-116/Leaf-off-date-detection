# Leaf-off-date-detection
Code of Leaf-off date detection
Four files are contained. 
_MODIS_dataprocess_new_ process the MOD13A1 data and result in high-quality EVI data.
_NTLdataprocess_ output the HDNTL data as daily file without missing days.
_leaf_off_examples_0904_ provides the general code for extracting FOD for a single pixel.
_leaf_off_cities_0904_ provides the general code for extracting FOD in city level.
The input dataset includes MOD13A1 (https://ladsweb.modaps.eosdis.nasa.gov/missions-and-measurements/products/MOD13A1)
and HDNTL (https://doi.org/10.5281/zenodo.14992989)
reference for methods
Chen, J., Jönsson, Per., Tamura, M., Gu, Z., Matsushita, B., & Eklundh, L. (2004). A simple method for reconstructing a high-quality NDVI time-series data set based on the savitzky–golay filter. Remote Sensing of Environment, 91(3–4), 332–344. https://doi.org/10.1016/j.rse.2004.03.014
Pei, Z., Zhu, X., Hu, Y., Chen, J., & Tan, X. (2025). A high-quality daily nighttime light (HDNTL) dataset for global 600+ cities (2012–2024). Copernicus GmbH. https://doi.org/10.5194/essd-2025-142
Zhang, X., Friedl, M. A., Schaaf, C. B., Strahler, A. H., Hodges, J. C. F., Gao, F., Reed, B. C., & Huete, A. (2003). Monitoring vegetation phenology using MODIS. Remote Sensing of Environment, 84(3), 471–475. https://doi.org/10.1016/s0034-4257(02)00135-9
