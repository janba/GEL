echo "Building libGEL.so"
python buildGEL.py
echo "Installing into /usr/local/lib"
sudo cp -p build/GEL/libGEL.so /usr/local/lib
cd ../src/GEL
sudo mkdir -p /usr/local/include/GEL
find . -name "*.h" -type f | xargs -I {} sudo cp --parents {} /usr/local/include/GEL
cd ../../GEL_UNIX


echo "Building libPyGEL.so"
python buildPyGEL.py
USER_SITE_DIR=$(python3 -m site --user-site)
echo "Installing into $USER_SITE_DIR"
mkdir -p ${USER_SITE_DIR}/PyGEL
cp -p build/PyGEL/libPyGEL.so ${USER_SITE_DIR}/PyGEL
cp -p ../src/PyGEL/*.py ${USER_SITE_DIR}/PyGEL
