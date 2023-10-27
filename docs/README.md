# Documentation reference guide

The RHM documentation is released at [github pages](https://silence2107.github.io/Rotochemical-heating-manager) 
and is powered by [jupyter-book](https://jupyterbook.org/en/stable/intro.html). In order to conduct changes to the documentation,
extra steps are required
  - Install jupyter-book, ghp-import pip packages
  - In this folder, run
```{bash}
# jupyter-book handles the process of creating webpage out of markdown files
jupyter-book build .
# ghp-import is a tool to push docs to a branch on GitHub and consecutively activate remarkup of Pages
ghp-import -n -p -f _build/html
```
