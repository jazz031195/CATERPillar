# CATERPillar

## **Introduction**

CATERPillar (**Computational Axonal Threading Engine for Realistic Proliferation**) is an advanced computational framework designed to simulate natural axonal growth using overlapping spheres as fundamental building blocks. By employing a biologically inspired approach, CATERPillar enables parallel axon development while effectively preventing collisions, allowing users to control key structural parameters such as **density, tortuosity, and beading**.

What sets CATERPillar apart is its ability to generate not only realistic **axonal architectures** but also **glial cell structures**, significantly enhancing the biological fidelity of tissue microstructure simulations. This makes it a valuable tool for studying brain tissue models and validating diffusion-based imaging techniques.

-----

## **How to Build and Run the Program**

CATERPillar uses CMake, making it fully cross-platform for Linux, macOS, and High-Performance Computing (HPC) environments.

### **1. Install Prerequisites**

Make sure you have a modern C++ compiler, CMake, Qt5, and Boost installed on your system.

**For Ubuntu / Debian Linux:**

```bash
sudo apt update
sudo apt install cmake qtbase5-dev libqt5datavisualization5-dev libboost-all-dev
```

**For macOS (using Homebrew):**

```bash
brew install cmake qt5 boost
```

### **2. Compile the Source Code**

Navigate to the downloaded repository folder in your terminal and run the following commands to generate the build files and compile the executable:

```bash
mkdir build
cd build
cmake ..
make
```

### **3. Run CATERPillar**

CATERPillar can be run in two distinct modes: an interactive graphical interface (ideal for exploring parameters) or a headless command-line mode (ideal for High-Performance Computing clusters and automated batch processing).

**Mode A: Graphical User Interface (GUI)**
Launch the GUI by running:

```bash
./CATERPillar
```
The GUI allows you to configure biophysical parameters and generate realistic numerical substrates. To generate a new substrate, click **"Grow Substrate"**. To visualize a previously generated substrate, click **"Visualise Substrate"**.

**Mode B: Headless / HPC Mode (JSON Configuration)**
To bypass the GUI entirely and run simulations in the terminal, pass the `--config` argument followed by the path to your JSON settings file:

```bash
./CATERPillar --config example_config.json
```
This mode allows you to execute heavy math and generate `.csv` files on remote servers without display capabilities (X11 forwarding is not required). An example configuration file (`example_config.json`) is included in the repository to help you correctly format your parameters.

-----

## **List of Parameters**

### **General Parameters:**

  * **Voxel Edge Length (μm):** Defines the length of the cubic substrate, determining the overall volume of the simulation space.
  * **Overlapping Factor:** Controls the spacing between consecutive spheres during axonal growth. The distance is computed as $\max(R_1, R_2) / F$, where $R_1$ and $R_2$ are the radii of two consecutive spheres, and $F$ is the chosen overlapping factor. A higher overlapping factor results in closer sphere placement, reducing the gaps. A recommended value for optimal results is **4**.
  * **Minimum Sphere Radius (μm):** Specifies the smallest allowable sphere radius within the voxel. This constraint helps prevent excessively narrow spaces that could impede Monte Carlo Simulations.
  * **Visualise Voxel:** Ticking this box allows the GUI to plot the 3D substrate once the growth is completed.

### **Axon Parameters:**

  * **Axons with Myelin ICVF (%):** Defines the volume fraction occupied by myelinated axons, including the myelin compartment. The g-ratio for each axon is assigned based on a log-linear relationship between the inner radius and the g-ratio.
  * **$K_1$, $K_2$, and $K_3$ parameters:** Defines the relationship between the inner diameter of the axons ($D_{\text{in}}$) and the myelin thickness:
    $$\text{Myelin thickness}=K_1+K_2\cdot D_{\text{in}}+K_3\cdot\log(D_{\text{in}})$$
    This log-linear fit was proposed by Lee et al. (DOI: 10.1007/s00429-019-01844-6) using $K_1=0.35$, $K_2=0.006$, and $K_3=0.024$.
  * **Axons ICVF (%):** Specifies the volume fraction of non-myelinated axons within the voxel.
  * **Number of Threads:** Determines the number of axons that can grow simultaneously during the simulation, impacting computational efficiency.
  * **Tortuosity ($\epsilon$):** Represents the standard deviation of the Gaussian distribution governing the 3D positioning of spheres during axonal growth. Higher values result in increased axonal tortuosity.
  * **Fibre Orientation Dispersion Function ($c_2$):** Defined as $\langle \cos^2\psi \rangle$, where $\psi$ is the angle between the axon growth direction and the z-axis. This parameter quantifies the degree of fibre orientation dispersion within the substrate.
  * **Number of Axon Populations:** Specifies the number of distinct axonal populations that can grow within the substrate (range: 1-3). Each population adopts a primary orientation perpendicular to the others. If two populations are selected, users can choose between a **sheet configuration** or an **interwoven configuration**.
  * **Beading Amplitude ($A$):** Defines the amplitude of axonal beading as a fraction of the axon's initial radius ($R$), influencing morphological variability. The radius of an axon changes with its length stochastically. The next sphere's radius is drawn from a normal distribution centered on the previous radius. If the computed radius falls outside the boundary $R \pm A \cdot R$, the normal distribution will instead be centered at that boundary value.
  * **Gamma Distribution for Radii ($\alpha$):** Shape parameter for the Gamma distribution governing axon radii. A recommended value for realistic axon widths is **4**.
  * **Gamma Distribution for Radii ($\beta$):** Scale parameter for the Gamma distribution governing axon radii. A recommended value for realistic axon widths is **0.25**.

### **Glial Cell Parameters:**

*(Note: Two different populations of glial cells can co-exist in a voxel. The GUI provides independent parameter boxes for each population.)*

  * **Somas ICVF (%):** Defines the volume fraction occupied by astrocyte somas within the substrate.
  * **Processes ICVF (%):** Specifies the volume fraction of astrocyte processes, contributing to the extracellular microstructure.
  * **Soma Radius Mean (μm):** The soma for each glial cell is drawn from a normal distribution with this value as the mean.
  * **Soma Radius Standard Deviation (μm):** The standard deviation used for the normal distribution of the soma radius.
  * **Mean Process Length (μm):** Represents the average length of processes within the substrate, measured from the soma to the tip of the process. Secondary processes (branching from primary ones) will have a shorter length, as the algorithm accounts for the distance already grown before the bifurcation.
  * **Standard Deviation for Process Length (μm):** Defines the variability in process lengths, contributing to the heterogeneity of the substrate.
  * **Number of Primary Processes:** The number of processes that emerge directly from the soma.
  * **Can glial cell population have branching:** If ticked, the glial cell can grow higher-order (secondary/tertiary) processes. If unticked, only primary processes will be present. If the target ICVF for processes cannot be reached with the selected number of primary processes, the system will automatically grow more of them.

-----

## **Output**

Once the simulation is complete, the **GUI** will display the generated substrate, allowing users to visually inspect the 3D microstructure. Additionally, the GUI provides analytical tools to plot:

  * **Axonal radii distribution**
  * **Axonal tortuosity distribution**
  * **Sholl analysis for astrocytes**

### **Generated Files**

Upon completion, the following output files are automatically saved in your selected directory:

1.  **`voxel.csv`** – Contains precise spatial and morphological data for the substrate, with each line defining a single sphere’s properties.
2.  **`voxel_info.txt`** – Provides metadata, summary statistics, and a record of the exact parameters used during the simulation.

The **`voxel.csv`** file is structured with the following columns for external analysis or Monte Carlo integration:

  * **`cell_type`:** The classification of the cell (`axon`, `glial_cell`, or `neuron`).
  * **`cell_id`:** The unique identifier of the cell (starting from 0).
  * **`component`:** The structural component (`soma`, `branch`, `axon`, or `spine`).
  * **`component_id`:** The identifier of the specific component (starting from 0).
  * **`X Y Z`:** The 3D spatial coordinates of the sphere center (μm).
  * **`inner_radius`:** The inner radius of the sphere (μm).
  * **`outer_radius`:** The outer radius of the sphere (μm). This will differ from the inner radius if myelin is present.

-----

## **Citation**

If you use CATERPillar in your research, please cite the following work:

> Nguyen-Duc JK, Brammerloh M, Cherchali M, De Riedmatten I, Perot JB, Rafael-Patino J, Jelescu IO. **CATERPillar: A Flexible Framework for Generating White Matter Numerical Substrates with incorporated Glial Cells.** *bioRxiv*, 2025.