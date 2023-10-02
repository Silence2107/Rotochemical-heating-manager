# jupyter-book handles the process of creating webpage out of markdown files
jupyter-book build .
# ghp-import is a tool to push a directory to a branch on GitHub
ghp-import -n -p -f _build/html