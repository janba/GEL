python buildGEL.py
sudo cp -p build/libGEL.a /usr/local/lib
cd ../src/GEL
sudo mkdir -p /usr/local/include/GEL
find . -name "*.h" -type f | xargs -I {} sudo cp --parents {} /usr/local/include/GEL
cd ../../GEL_UNIX

python buildPyGEL.py

USER_SITE_DIR=$(python3 -m site --user-site)
mkdir -p ${USER_SITE_DIR}/PyGEL
cp -p build/libPyGEL.so ${USER_SITE_DIR}/PyGEL
cp -p ../src/PyGEL/*.py ${USER_SITE_DIR}/PyGEL
