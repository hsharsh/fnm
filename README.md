The program has been written and tested with C++17
One can compile the files using the following syntax:

g++ -std=c++17 (input-file) -o (output-file)

Eigen is used for major linear alegbra. Version 3.3.9 is need for some implementations. So download the same from the following page. Change the name of the folder in the include statement based on the folder you place the files in the /usr/include
http://eigen.tuxfamily.org/index.php?title=Main_Page

Note that while using operations like x = x(slice), it is possible to get errors because the () operator uses a reference to the variable and thus the above operation ends up reading and writing from the same variable causing errors.

Reading csv files is taken from the following answer on Stackoverflow:
https://stackoverflow.com/a/39146048

Reading from a configuration file taken from the folowing blogpost:
https://www.walletfox.com/course/parseconfigfile.php