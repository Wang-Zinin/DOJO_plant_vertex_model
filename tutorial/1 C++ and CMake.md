# How to Use C++ and CMake
This tutorial will tell you about how to use C++ and CMake, which are prerequisites for the PlantVertexModel simulation. For those who are already familiar with C++ and CMake, please skip this tutorial.

This tutorial is a very simple one and does not aim to tell the details of C++ programming. If you want to learn more about C++ coding, please read "C++ Primer Plus" (ISBN-13: 978-0-321-77640-2).

# For Windows Users 
There are several ways to compile and build C++ files and use CMake. We will use the WSL (Windows Subsystem for Linux) ubuntu here.  
## Installing WSL
Open PowerShell or Windows Command Prompt (you may find them through search "PowerShell"/"Windows Command Prompt" in your taskbar) in adminstrator mode by right-clicking and selecting "Run as administrator", enter: 
```PowerShell
wsl --install
```
This command will enable the features necessary to run WSL and install the Ubuntu distribution of Linux. Then please restart the computer. 

The first time you open ubuntu, it will require you to register a user's name and password. Please do the register and remember your password. Installing or updating packages sometimes require the password.

To check if you have successfully installed WSL, please enter the following command in PowerShell or Windows Command Prompt:
```PowerShell
wsl -l -v
```
This command will show you the existing wsl version you have.

For details of wsl please refer to  https://learn.microsoft.com/en-us/windows/wsl/install

## Install Visual Studio Code (optional)

Visual Studio Code is a good code editor. You can use it to modify codes. You can install it directly from Microsoft store. Then you can open the codes files or file folders by Visual Studio Code.  

## Install C++ compiler (g++) and CMake in WSL Ubuntu

Please open the WSL Ubuntu and enter the following codes 

```Ubuntu
sudo apt-get update
sudo apt install build-essential -y
sudo apt-get install cmake
```
The C++ compiler is in "build-essential". 

To check if you have successfully installed C++ compiler (g++) and CMake. Please enter the following command in ubuntu:

```Ubuntu
homebrew 
g++ --version
cmake --version
```
These commands will show you the versions of your current g++ and CMake. 

# For Mac Users
## Install HomeBrew
Homebrew is a package manager for macOS. It can install the stuff (packages) you need. We will use HomeBrew to install C++ compiler and CMake in your macOS.  
Open you macOS terminal, then enter the following command: 
```macOS
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

To check if you have successfully installed HomeBrew, please enter the following command: 
```macOS
brew --version
```

Reference:  https://brew.sh/

## Install C++ compiler and CMake by HomeBrew 
To install C++ compiler and CMake by HomeBrew, please enter the following commands in macOS terminal: 

```macOS
homebrew install gcc
homebrew install cmake
```

To check if you have successfully installed gcc and cmake, please enter the following command: 

```macOS
gcc --version
cmake --version
```

# Hello World in C++
Open your ubuntu/macOS terminal, go to the desired directory by using the following commands:

```
ls                #show the list of files and file directory in the current working directory
cd FileDirectory  #Navigate to the desired directory
```

Then create a file named "HelloWorld.cpp" in the desired directory. This could be easily done by Visual Studio Code or though the command in ubuntu/macOS terminal:
```
touch HelloWorld.cpp #create a file named "HelloWorld.cpp"
```

Then open the HelloWorld.cpp file by either Visual Studio Code or by:
in ubuntu:
```
nano HelloWorld.cpp
```
in mac:
```
open HelloWorld.cpp
```

Then enter the following codes in the "HelloWorld.cpp" file:  
```C++
#include <iostream>                             //a Preprocessor directive. 
                                                //The standard library named <iostream> is used. 
                                                //This library includes functions to input and output. "io" means input and output. 
                                                //In this program, we will cout to print the HelloWorld message on our terminal.

int main()                                      //function header
{                                               //start of function body
    std::cout<<"Hello World C++ !"<<std::endl;  //std:: means standard library; cout is a function to print message on terminal
    return 0;                                   //terminate main()
}                                               //end of function body
```

Then in ubuntu/macOS terminal, enter command: 
```
g++ HelloWorld.cpp
./a.out
```
There will be a "Hello World C++ !" message printed out in your terminal.

Congradulations ! You have compiled a simple C++ code.

Here is a flowchart to explain what g++ compiler have done: 
<img src="../flowchart_helloworld.jpg" width="800"/>

# Hello World in CMake

Create a file named "CMakeLists.txt" in the same file folder of "HelloWorld.cpp" file and enter the following codes:
``` CMake
cmake_minimum_required(VERSION 3.10)
project(hello)
set(CMAKE_CXX_STANDARD 11)
add_executable(
    ${PROJECT_NAME}
    HelloWorld.cpp
)
```
Then in ubuntu/macOS terminal, enter command: 
```
mkdir build
cd build
cmake ..
make
./hello
```
There will be a "Hello World C++ !" message printed out in your terminal.

Congradulations ! You have compiled C++ files through CMake. 
