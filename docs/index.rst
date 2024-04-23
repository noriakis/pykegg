pykegg
===================================

**pykegg** is a Python library for visualization and analysis of KEGG information.


Motivation
----------

KEGG information (especially KGML) can be parsed effortlessly by the great package `biopython`. The package has a built-in module `KGML_vis <https://biopython.org/docs/1.76/api/Bio.Graphics.KGML_vis.html>`_. Also, there are great tutorials for visualizing KEGG information in Python, like `this notebook <https://nbviewer.org/github/widdowquinn/notebooks/blob/master/Biopython_KGML_intro.ipynb>`_ by Dr. Leighton Pritchard, and `this tutorial <https://ctskennerton.github.io/2020/09/27/drawing-kegg-pathway-maps-using-biopython-and-matplotlib/>`_ by Dr. Connor Skennerton. There is a nice Python library for visualization of KEGG information with enrichment analysis, `keggtools <https://keggtools.org/>`_. Extending these packages, I would like to create a package that can address the following problems.

* Parse the information from KEGG as a network and conduct network analysis.
* Combine with the other omics analysis libraries like ``scanpy``, ``PyDESeq2``.
* Visualize the KEGG information by libraries like ``OpenCV`` or ``plotnine`` in easily customizable ways, using grammer of graphics.

For R environment, please refer to `ggkegg <https://noriakis.github.io/software/ggkegg>`_.

Analysis
--------
Please see :ref:`notebooks <notebooks>` for the visualization example and the example analyses using ``PyDESeq2``, ``scanpy``, and ``GSEApy``.


.. note::

   This project is under active development.

.. warning::

    **pykegg** uses KEGG API, and is restricted to academic use by academic users belonging to academic institutions.

.. toctree::
  :caption: Getting started
  :hidden:

  usage
  modules

.. toctree::
  :caption: Example analysis
  :hidden:

  example_analysis
