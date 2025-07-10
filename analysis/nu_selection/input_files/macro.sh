g++ -std=c++17 -o generate_input_list generate_input_list.cc
if [ $? -eq 0 ]; then
    echo "Compilation successful"
else
    echo "Compilation failed. "
fi