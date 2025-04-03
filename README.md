# 🚀 Project Name: R Environment Setup with `renv`

This project uses [`renv`](https://rstudio.github.io/renv/) for package management to ensure **consistent environments** across different systems (Windows & macOS).

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
Follow these steps to get started after cloning the repository.

### 1️⃣ Install `renv` (if not installed)
Open **R/RStudio** and run:  
```r
install.packages("renv")  # Installs renv package
```

### 2️⃣ Restore the Project Environment
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

### 🖥 Windows Users
- If package installation fails, make sure you have **Rtools** installed:  
  Download from [CRAN Rtools](https://cran.r-project.org/bin/windows/Rtools/)  

### 🍏 Mac Users
- If Bioconductor packages fail, try:  
  ```r
  BiocManager::install("your_package_name", force = TRUE)
  ```

---

## 📂 Files Included in This Repo

| File | Description |
|------|------------|
| `renv.lock` | Stores package versions for reproducibility |
| `renv/activate.R` | Ensures `renv` is activated in the project |
| `.gitignore` | Ignores unnecessary files like `renv/library/` |

---

## 📢 Need Help?
If you run into issues, feel free to ask! 🚀  
