PYTHONWARNINGS="ignore"
rm -r .local
python -m pip install --no-cache-dir --user --upgrade pip  2>/dev/null
python -m pip install --no-cache-dir --user python-dateutil==2.8.2  2>/dev/null
python -m pip install --no-cache-dir --user basemap==1.3.9  2>/dev/null
python -m pip install --no-cache-dir --user shapely==1.8.5 --no-binary shapely --ignore-installed  2>/dev/null
python -m pip install --no-cache-dir --user wrf-python  2>/dev/null
python -m pip install --no-cache-dir --user metpy  2>/dev/null
python -m pip install --no-cache-dir --user cartopy  2>/dev/null
python -m pip install --no-cache-dir --user netCDF4  2>/dev/null
#setfacl -R -m d:u:barte035:rwx .
#setfacl -R -m u:barte035:rwx .
