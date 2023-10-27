# Documentation reference guide

The RHM documentation is released at [github pages](https://silence2107.github.io/Rotochemical-heating-manager) 
and is powered by [jupyter-book](https://jupyterbook.org/en/stable/intro.html). In order to conduct changes to the documentation,
extra steps are required
  - Install jupyter-book, gh-pages pip packages
  - In this folder, run
```{bash}
# jupyter-book handles the process of creating webpage out of markdown files
jupyter-book build .
# ghp-import is a tool to push a directory to a branch on GitHub
ghp-import -n -p -f _build/html
```
