# TOG Workshop Apr 15

Core package dependencies:

- tidyverse
- patchwork
- UpSetR
- bioconductor
- ExperimentHub
- limma
- muscData
- muscat
- scater


## 📥 Setup Instructions
Follow these steps to get started.

### 1️⃣ Clone the Repository
You can clone the repository using either **RStudio** or the command line.

#### Option 1: Using RStudio
1. Open **RStudio**.
2. Go to **File** → **New Project** → **Version Control** → **Git**.
3. In the **Repository URL** field, enter:
   ```sh
   https://github.com/tonyliang19/tog-workshop
   ```
4. Choose a local directory to store the project.
5. Click **Create Project** to clone the repository.

#### Option 2: Using the Command Line
Alternatively, you can clone the repository with:
```sh
git clone https://github.com/tonyliang19/tog-workshop
```

Then, open **R/RStudio** and navigate to the project directory.

---

### 2️⃣ Install `renv` (if not installed)
Open **RStudio** and run:
```r
install.packages("renv")  # Installs renv package
```

### 3️⃣ Restore the Project Environment
Run this in R to install all required packages:
```r
renv::restore()  # Installs all packages from renv.lock
```
This ensures you have the **exact package versions** used in this project.

---

## 🔧 Additional Setup for macOS Users
Some Bioconductor packages require system dependencies. If you are on macOS, install them using **Homebrew**:

```sh
brew install gfortran
brew install libxml2
```
If you don’t have **Homebrew**, install it first:  
```sh
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

---

## 🐞 Troubleshooting

- You need R version of at least 4.3 and Bioconductor >= 3.18


### 🖥 Windows Users

- If package installation fails, make sure you have **Rtools** installed:  
  Download from [CRAN Rtools](https://cran.r-project.org/bin/windows/Rtools/)

### 🍏 Mac Users
- If Bioconductor packages fail, try:  
  ```r
  BiocManager::install("your_package_name", force = TRUE)
  ```

## 📢 Need Help?
If you run into issues, feel free to ask! 🚀

