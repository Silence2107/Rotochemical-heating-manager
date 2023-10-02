# First run

RHM follows standard _bin/include/src_ structure, with program's logic assembled under _project/_ folder. _3rd-party/_ folder contains several IO libraries (ready-to-use), and _presupplied/_ contains presupplied data. _Makefile_ is also supplied for compilation.

Since there are many different main programs under _project/_, some of which may run with settings provided by the user and some may not, making executables is yet user's responsibility.

To perform a test run with presupplied data, run make on some typical script, e.g.

```bash
make release app=project/M-R_diagram/m_r_diagram
```

<!-- If successful, the corresponding binary will be put under \textit{bin/}, following same path as it took to the app. In this case, executing binary would work like follows

\begin{lstlisting}
    bin/project/M-R_diagram/m_r_diagram.out --help
\end{lstlisting}

\textit{$--$help} invokes manual message for all standardized RHM programs.

To finally see whether physics is in order on your machine, this binary (at least the way it is provided) must be supplied with inputfile. 

\textbf{Q: }Given your knowledge from binary's manual, supply it with inputfile, which lies under \textit{saved/essentials/RHMconfig.json}.

\textbf{Expected output :} 
\begin{itemize}
    \item $(M, R)$ pairs of order $(2 \text{M}_\odot, 10 \text{km})$ based on various center density fractions.
    \begin{itemize}
        \item phenomenal work
    \end{itemize}
\end{itemize}

\textbf{Unexpected output : }
\begin{itemize}
    \item $(M, R)$ pairs of $(0 \text{M}_\odot, 0 \text{km})$ \begin{itemize}
        \item inputfile is not supplied
    \end{itemize}
    \item "(..) Cannot open file (..)" 
    \begin{itemize}
        \item inputfile path is supplied, but is not recognized as valid
        \item Check path's spelling against \textit{saved/essentials/RHMconfig.json}
        \item If correct, make sure there exists a valid file under \\ \textit{saved/EoS\_bank/APR\_EOS\_Acc\_Fe\_RHMstandard.dat}
        \item If exists, make sure there's an entry in \\
        \textit{saved/essentials/RHMconfig.json} ["EoSSetup"]["Datafile"]["Path"] that leads to the file above. Specify absolute path if in doubts.
    \end{itemize}
    \item "keyword argument (..) must have value."
    \begin{itemize}
        \item program's key is supplied, but actual value to it is not
    \end{itemize}
    \item Whatever else happened
    \begin{itemize}
        \item Tell me I messed
    \end{itemize}
\end{itemize}
\textbf{A: }
\begin{lstlisting}
    bin/project/M-R_diagram/m_r_diagram.out --inputfile saved/essentials/RHMconfig.json
\end{lstlisting} -->

If successful, the corresponding binary will be put under _bin/_, following same path as it took to the app. In this case, executing binary would work like follows

```bash
bin/project/M-R_diagram/m_r_diagram.out --help
```

`--help` invokes manual message for all standardized RHM programs.

To finally see whether physics is in order on your machine, this binary (at least the way it is provided) must be supplied with inputfile.

**Q:** Given your knowledge from binary's manual, supply it with inputfile, which lies under _presupplied/Inputfile/RHMconfig.json_.

**Expected output :**
- $(M, R)$ pairs of order $(2 \text{M}_\odot, 10 \text{km})$ based on various center density fractions.
    - phenomenal work

**Unexpected output :**
- $(M, R)$ pairs of $(0 \text{M}_\odot, 0 \text{km})$ 
    - inputfile is not supplied
- "(..) Cannot open file (..)"
    - inputfile path is supplied, but is not recognized as valid
    - Check path's spelling against _presupplied/Inputfile/RHMconfig.json_
    - If correct, make sure there exists a valid file under _saved/EoS\_bank/APR\_EOS\_Acc\_Fe\_RHMstandard.dat_
    - If exists, make sure there's an entry in _presupplied/Inputfile/RHMconfig.json_ `[\"EoSSetup\"][\"Datafile\"][\"Path\"]` that leads to the file above. Specify absolute path if in doubts.
- "keyword argument (..) must have value."
    - program's key is supplied, but actual value to it is not
- Whatever else happened
    - Tell me I messed

**A:**
```bash
bin/project/M-R_diagram/m_r_diagram.out --inputfile presupplied/Inputfile/RHMconfig.json
```