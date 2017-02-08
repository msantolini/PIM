#PIM 

The Pairwise Interaction Model (PIM) describes the binding of a Transcription
Factor (TF) to its DNA binding sites (BS) using the statistics of both single
nulceotides and pairs of nucleotides. In the physics jargon, this model is
known as an inhomogeneous Potts models with 4 colors.

##DESCRIPTION


The `src` folder contains scripts to estimate a variety of binding models from
ChIPseq or ChIP-chip data, namely:

* The **Position Weight Matrix (PWM) model**, where the TF binds to each of the
  nucleotides independently.

* The **PWM mixture model**, where the TF is assumed to adopt several binding
  conformations, each of which is represented by a PWM and has a certain
  probability to occur.

* The **Pairwise Interactions Model (PIM)**, where pairwise contributions from
  pairs of nucleotides are taken into account.

* The **Nearest Neighbor Model (NNM)**, which is a PIM with pairwise
  contributions restricted to nearest-neighbor nucleotides.

The `data` folder contains an example ChIP-chip dataset for the Drosophila TF
Twist that is used as default for the models estimation. It contains two
folders:

* The `data/nmers-chips` folder contains a file giving the N-mers
(here N=12) to be searched for binding sites. For each ChIP peak, a unique
identifier was assigned (first column) and all N-mers from both strands were
considered (second column). When learning a model, only the 50% best peaks
having one or several high scoring BS are taken into account. The identifier
allows to keep track of the peak on which the BS resides. Note that nothing
restricts the models to be learnt from ChIP-chip or ChIPseq data, and such
an N-mers file could be generated from other kinds of experiments.

* The `data/motinits` folder contains a file giving an initial frequency matrix
  for the TF of interest, where rows correspond to A,T,C,G frequencies and
  columns to different positions of the site.

For a given TF, both files should have the same name, which correspond the the
`MAT` variable in the matlab scripts. For example, in the default case,
`MAT='twi'` and the files are named `twi.dat`.


##FIRST STEPS

Most of the code is written in matlab. It has been executed with matlab
versions greater than R2011b but has not been tested with older versions.

Before using the code, two things need to be taken care of.

1. First, the exhaust.c file has to be compiled using the [mex
   function](http://www.mathworks.com/help/matlab/ref/mex.html). To do that,
   open matlab in the `src` directory and execute

    ```matlab
    mex exhaust.c
    ```

    This will produce a matlab compatible file that lies at the heart of the PIM
    estimation. The exhaust function loops across al possible N-mers (in our case,
    N=12) to compute several quantities of interest, such as the single and
    pairwise nucleotides probabilities estimated by the model or the Partition
    Function used for probability normalization. The computation time therefore
    grows exponentially with N, which prohibits from using large values (N=12
    has proven to be a good compromise).

2. Second, the runMefirst.sh script should be executed:

    ```sh
    ./runMeFirst.sh
    ```

    This will produce the `src/load_paths.m` file that defines the **absolute** paths
    of the `data` and `src` folders, as well as the `results` directory that will
    be created when running the matlab scripts. These paths can be adapted to your
    need by modifying the definitions of the SRCDIR, DATADIR and RESDIR
    variables in the runMeFirst.sh script. All these directories are located at
    the root of the PIM folder per default. **If the name or the location of
    any of the parent folders changes, this script should be executed again.**

##USAGE

The `src/compAllModels.m` script can be executed to learn the different models on the
sample Twist data. To do so, open matlab in the `src` folder and execute:

```
compAllModels
```

The different models can also be estimated separately. For example, if one is solely interested in learning the PWM model, then the `src/PWM_model.m` script should be executed.

The important parameters, such as the length of the binding sites N or the
motif name, are all located in the file `src/load_params.m`, which is loaded in
the beginning of all scripts. This is were any change should be made if one
wants to conduct a custom study (look at a different TF, change N, etc).


##CITATION

Santolini M, Mora T, Hakim V (2014) A General Pairwise Interaction Model Provides an Accurate Description of In Vivo Transcription Factor Binding Sites. PLoS ONE 9(6): e99015. doi:10.1371/journal.pone.0099015
