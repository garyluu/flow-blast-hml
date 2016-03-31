Installation
=============


Nextflow
------------------------

.. code-block:: shell

	curl -fsSL get.nextflow.io | bash


Docker
------------------------

.. code-block:: shell

	docker-machine start default
	eval $(docker-machine env default)


Running the tutorial
------------------------

.. code-block:: shell

	git clone https://github.com/nmdp-bioinformatics/flow-blast-hml
	cd flow-blast-hml
	./nextflow run nmdp-bioinformatics/flow-blast-hml -with-docker \
	nmdpbioinformatics/docker-blast-hml \
	--hml tutorial/ex00_ngsp_expected.xml \
	--outdir tutorial --report 0
