echo "Building suk"
cd suk
mkdir -p build
cd build
rm -r *
cmake ..
make -j 10
make
cd ../..

echo "Building overlap"
cd overlap
./install_deps.sh
mkdir -p build
cd build
rm -r *
cmake ..
make -j 10
make
cd ..

echo "Building hypo polisher"
cd polisher
./install_deps.sh
mkdir -p build
cd build
rm -r *
cmake ..
make -j 10
make
cd ../..

echo "Building scaffolder"
cd scaffold
./install_deps.sh
mkdir -p build
cd build
rm -r *
cmake ..
make -j 10
make
cd ../..

mkdir -p run_all
cp run_all.sh run_all/
cp suk/build/bin/suk run_all/
cp overlap/run_overlap.sh run_all/
cp overlap/join_overlap.py run_all/
cp overlap/filter_overlap.py run_all/
cp polisher/build/bin/hypo run_all/
cp scaffold/run_overlap.sh run_all/
cp scaffold/join_overlap.py run_all/
cp scaffold/filter_overlap.py run_all/
