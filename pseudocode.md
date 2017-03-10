~~Strike through~~
[] task lists

```javascript
function
```

	code1
	code2
	function

%=========================================================
% Titel

\title{SSimp user guide}
\author{ Aaron McDaid\thanks{\texttt{<aaron.mcdaid@gmail.com>}}
 \and Zolt{\'a}n Kutalik\thanks{\texttt{<zoltan.kutalik@unil.ch>}}
\and Sina R\"{u}eger\thanks{\texttt{<sina.rueeger@gmail.com>}}}
\date{\today\\for \ssimp{} version \qversion}


%=========================================================

\begin{document}

\maketitle

\clearpage
\tableofcontents

%=========================================================
% installing ssimp
\clearpage
\section{Installing or compiling SSimp}



%=========================================================
% Simple default usage
\section{Quick start}

\begin{verbatim}
ssimp --data gwas.txt --out qres.txt
\end{verbatim}

This just imputes all SNPs that are within the gwas SNPs and in the reference panel.

\begin{verbatim}
ssimp --data gwas.txt --snp file.with.all.snps.to.impute.txt --out qres.txt
\end{verbatim}
This will only impute SNPs mentioned in snp.



%=========================================================
% Parameters
\section{Arguments}

\begin{description}
\item [\dash{data}] quicktest or snptest format or own format: identifier(SNP or Pos:Chr on hg18),
  Z, N, A1, A2. Z= Z statistics, sample size N. If sample size is NA, Z will be imputed. Otherwise the standardized
  effect size. If P, beta are provided, Z will be computed. If Z nor N is available, bst is
  searched. A1 (reference allele) A2 (effect allele) can be in upper or lower case.
\item [\dash{lambda}] sdfsdf
\item [\dash{missings}] truncate SNPs
\item [\dash{out}] sdfsdf
\item [\dash{snps}] file with all snps to impute. no header, will guess wether SNP or SNP2. 
\item [\dash{ref}] reference panel, either uk10k, 1kg, or own (see further down)
\item [\dash{ref.pop}] subpopulation (in case of 1kg) [tsi, lkkjlk]
\item [\dash{ref.allelematching}] allele matching
\end{description}

\subsection{ref}
Provide your own vcf data, and the path to it.

\subsection{identifier}
working with SNP2

%=========================================================
% 
\section{Output}
There is an \texttt{.log} file providing the output and possible warning messages. The .out file has
the following columns:

\begin{description}
\item [\dash{SNP}] SNP
\item [\dash{Chr}] Chromosome (only 1 to 22 right now)
\item [\dash{Pos}] Position (HG18?)
\item [\dash{SNP2}] Chr:Pos
\item [\dash{bst}] bst
\item [\dash{N.est}] estimation of N
\item [\dash{r2.pred}] impqual
\item [\dash{A1}] reference allele
\item [\dash{A2}] effect allele

\end{description}

\texttt{.out.not} lists all the SNPs that were not found in the reference panel 
\texttt{.out.tnot} lists all the tSNPs that were not used as tSNPs (but needed): not found in the reference panel or had an NA in Z

If SNP has an NA as bst, it means that there was no tagSNP present.

%=========================================================
% Anhang
\newpage
\appendix


%=========================================================
% Misc
\section{Edge-cases}
\label{misc}

\begin{itemize}
\item tSNP not in reference panel
\item SNP not in reference panel
\item ref alleles not swaped as in reference panel
\item A1 and A2 are not matching with reference panel
\item NA in Z of tSNP
\item two tSNPs or SNPs same Pos and chr
\item !is.na(SNP) but is.na(SNP2)
\item types of missings: NA, "-"

\end{itemize}

%--------------------------------------
\end{document}
%--------------------------------------
