mkdir -p output_files

# Compile from the 'sources' directory
g++ -std=c++17 -O2 -g -pedantic -Wall \
    sources/main.cc sources/utils.cc sources/cuts.cc sources/selection_input.cc \
    -o select_n \
    -Iheaders \
    `root-config --cflags --glibs` -lEG \
    -I/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_06_01b/include \
    -I$SRPROXY_INC \
    -L$DUNEANAOBJ_LIB -lduneanaobj_StandardRecordProxy \
    -lduneanaobj_StandardRecord \
    -lduneanaobj_StandardRecord_dict

if [ $? -eq 0 ]; then
    echo "Compilation successful"
else
    echo "Compilation failed"
    return 1 2>/dev/null || exit 1
fi
