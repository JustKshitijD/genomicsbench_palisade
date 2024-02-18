
.. _installation:

Getting Started
===============


Installation Methods
--------------------

Medaka can be installed in one of several ways.

**Installation with conda**

Perhaps the simplest way to start using medaka on both Linux and MacOS is
through conda; medaka is available via the
`bioconda <https://anaconda.org/bioconda/medaka>`_ channel:

.. code-block:: bash

    conda create -n medaka -c conda-forge -c bioconda medaka


**Installation with pip**

Medaka can be installed using the python package manager, pip:

.. code-block:: bash

    pip install medaka

On Linux platforms this will install a precompiled binary, on MacOS (and other)
platforms this will fetch and compile a source distribution.

We recommend using medaka within a virtual environment, viz.:

.. code-block:: bash

    virtualenv medaka --python=python3 --prompt "(medaka) "
    . medaka/bin/activate
    pip install medaka

.. note::

    Using this method requires the user to provide several binaries:

    `samtools <https://github.com/samtools/samtools>`_,
    `minimap2 <https://github.com/lh3/minimap2>`_,
    `tabix <https://github.com/samtools/htslib>`_,
    `bgzip <https://github.com/samtools/htslib>`_, and
    `bcftools <https://github.com/samtools/bcftools>`_,

    and place these within the `PATH`. `samtools` version 1.9, `bcftools` version 1.10.2 and `minimap2`
    version 2.17 are recommended as these are those used in development of
    medaka.


**Installation from source**

Medaka can be installed from its source quite easily on most systems.

.. note::

    Before installing medaka it may be required to install some
    prerequisite libraries, best installed by a package manager. On Ubuntu
    theses are:

    bzip2 g++ zlib1g-dev libbz2-dev liblzma-dev libffi-dev libncurses5-dev
    libcurl4-gnutls-dev libssl-dev curl make cmake wget python3-all-dev python-virtualenv

A Makefile is provided to fetch, compile and install all direct dependencies
into a python virtual environment. To setup the environment run:

.. code-block:: bash

    # Note: certain files are stored in git-lfs, https://git-lfs.github.com/,
    #       which must therefore be installed first.
    git clone https://github.com/nanoporetech/medaka.git
    cd medaka
    make install
    . ./venv/bin/activate

Using this method both ``samtools`` and ``minimap2`` are built from source and need
not be provided by the user.


**Using a GPU**

Since version 1.1.0 ``medaka`` uses Tensorflow 2.2, prior versions used Tensorflow 1.4.
For ``medaka`` 1.1.0 and higher installation from source or using ``pip`` can make
immediate use of GPUs. However, note that the ``tensorflow`` package is compiled against
specific versions of the NVIDIA CUDA and cuDNN libraries; users are directed to the
`tensorflow installation <https://www.tensorflow.org/install/gpu>`_) pages
for further information. cuDNN can be obtained from the
`cuDNN Archive <https://developer.nvidia.com/rdp/cudnn-archive>`_, whilst CUDA
from the `CUDA Toolkit Archive <https://developer.nvidia.com/cuda-toolkit-archive>`.

For ``medaka`` prior to version 1.1.0, to enable the use of GPU resource it is
necessary to install the ``tensorflow-gpu`` package. Using the source code from github
a working GPU-powered ``medaka`` can be configured with:


.. code-block:: bash

    git clone https://github.com/nanoporetech/medaka.git
    cd medaka
    sed -i 's/tensorflow/tensorflow-gpu/' requirements.txt
    make install

*GPU usage notes:*

Depending on your GPU, ``medaka`` may show out of memory errors when running.
To avoid these the inference batch size can be reduced from the default
value by setting the ``-b`` option when running ``medaka_consensus``. A value
``-b 100`` is suitable for 11Gb GPUs.

For users with RTX series GPUs it may be required to additionally set an
environment variable to have ``medaka`` run without failure:

.. code-block:: bash

    export TF_FORCE_GPU_ALLOW_GROWTH=true

In this situation a further reduction in batch size may be required.


.. _sequence_correction:

Sequence correction
-------------------

After installing the software (see :ref:`installation`), `medaka` can be run
using its default settings through the `medaka_consensus` program. An
assembly in `.fasta` format and basecalls in `.fasta` or `.fastq` format are
required (see :ref:`basecalling_and_draft_assembly` for an detailed example
of one method of obtaining these). More details and background can be found in
:ref:`draftorigin`.


The program uses both `samtools` and `minimap2`.
If medaka has been installed using the from-source method these will be present
within the medaka environment, else they will need to be provided by the user.

.. code-block:: bash

    source ${MEDAKA}  # i.e. medaka/venv/bin/activate
    NPROC=$(nproc)
    BASECALLS=basecalls.fa
    DRAFT=draft_assm/assm_final.fa
    OUTDIR=medaka_consensus
    medaka_consensus -i ${BASECALLS} -d ${DRAFT} -o ${OUTDIR} -t ${NPROC} -m r941_min_fast_g330

The variables ``BASECALLS``, ``DRAFT``, and ``OUTDIR`` in the above should be set
appropriately. When ``medaka_consensus`` has finished running, the consensus
will be saved to ``${OUTDIR}/consensus.fasta``.

.. warning::

  For best results it is important to specify the correct model, ``-m`` in the
  above, according to the basecaller used. Allowed values can be found by
  running ``medaka tools list_models``.


Medaka models are named to indicate i) the pore type, ii) the sequencing
device (MinION or PromethION), iii) the basecaller variant, and iv) the
basecaller version:

.. code-block:: bash

      {pore}_{device}_{caller variant}_{caller version}

For example the model named ``r941_min_fast_g303`` should be used with data from
MinION (or GridION) R9.4.1 flowcells using the fast Guppy basecaller version
3.0.3. By contrast the model ``r941_prom_hac_g303`` should be used with PromethION
data and the high accuracy basecaller (termed "hac" in Guppy configuration
files). Where a version of Guppy has been used without an exactly corresponding
medaka model, the medaka model with the highest version equal to or less than
the guppy version should be selected.


Improving parallelism
~~~~~~~~~~~~~~~~~~~~~

The ``medaka_consensus`` program is good for simple datasets but perhaps not
optimal for running large datasets at scale. examples. A higher level of
parallelism can be achieved by running independently the component steps
of ``medaka_consensus``. The program performs three tasks:

1. alignment or reads to input assembly (via ``mini_align`` which is a thin
   veil over ``minimap2``)
2. running of consensus algorithm across assembly regions
   (``medaka consensus``, note no underscore!)
3. aggregation of the results of 2. to create consensus sequences
   (``medaka stitch``)

The three steps are discrete, and can be split apart an run independently. In
most cases, Step 2. is the bottleneck and can be trivially parallelized. The
``medaka consensus program`` can be supplied a ``--regions``
argument which will restrict its action to particular assembly sequences from
the ``.bam`` file output in Step 1. Therefore individual jobs can be run for batches
of assembly sequences simultaneously. In the final step, ``medaka stitch``
can take as input one or more of the ``.hdf`` files output by Step 2.

So in summary something like this is possible:

.. code-block:: bash

    # align reads to assembly
    mini_align -i basecalls.fasta -r assembly.fasta -P -m \
        -p calls_to_draft.bam -t <threads>
    # run lots of jobs like this, change model as appropriate
    mkdir results
    medaka consensus calls_to_draft.bam results/contigs1-4.hdf \
        --model r941_min_fast_g303 --batch 200 --threads 8 \
        --region contig1 contig2 contig3 contig4
    ...
    # wait for jobs, then collate results
    medaka stitch results/*.hdf polished.assembly.fasta

It is not recommended to specify a value of ``--threads`` greater than 8 for
``medaka consensus`` since the compute scaling efficiency is poor beyond this.
Note also that ``medaka consensus`` may been seen to use resource equivalent to
``<threads> + 4`` as an additional 4 threads are used for reading and preparing
input data.

Human variant calling
---------------------

Variant calling on human samples can be performed with the `medaka_variant` workflow:

    medaka_variant -f <REFERENCE.fasta> -b <reads.bam>

which requires simply a path to the human reference sequence, and a `.bam` file
of reads aligned to the reference.

Bacterial (ploidy-1) variant calling
------------------------------------

Variant calling for monoploid samples is enabled through the `medaka_haploid_variant`
workflow:

    medaka_haploid_variant <reads.fastq> <ref.fasta>
    
which requires the reads as a `.fasta` or `.fastq` and a reference sequence as a
`.fasta` file.
