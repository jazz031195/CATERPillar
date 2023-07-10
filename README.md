# Sim_Growth
## Installations
### MacOs
Type the following in a terminal :
- /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
- brew install freeglut

## Parameters
The parameters are all defined in the main.cpp.
- number_axons (unsigned int): Number of axons to grow in the substrate.
- alpha (double) : alpha parameter for a Gamma Distribution. The radii of the axons arew taken from this distribution.
- beta (double) : beta parameter for a Gamma Distribution. The radii of the axons arew taken from this distribution.
- min_l (Eigen::Vector3d) : The coordinates of a "minimum point" in space to define the corner of the 3D box to work in.
- max_l (Eigen::Vector3d) : The coordinates of a "maximum point" in space to define the corner of the 3D box to work in.
- min_radius (double) : minimum radius of the axons.
- tortuous (bool) : if true, the axons will be tortuouss, overwise they will grow straight.

# Run the program
- In the compiler.sh, modify the following line : cd "/home/localadmin/Documents/Sim_Growth/src/", to the path that leads to the src folder on your local machine.
- If you are using a MacOs, change in compiler.sh the command line to : command="g++ -I/usr/include -O3 -std=c++11 -std=c++0x -o Growth-Animation.exe $src_files -I. -L/usr/local/lib -lsfml-graphics -lsfml-window -lsfml-system -lGL -lglut -lGLU"
- In the terminal, use "cd" to go to the Sim_Growth folder.
- Type : "./run.sh"
