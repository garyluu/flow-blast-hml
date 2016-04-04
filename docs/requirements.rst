Running the Tutorial
================================

Clone github repository
-----------------------
.. code-block:: shell

	git clone https://github.com/nmdp-bioinformatics/flow-blast-hml
	cd flow-blast-hml

Install nextflow
----------------
.. code-block:: shell

	curl -fsSL get.nextflow.io | bash

Start docker machine
--------------------
.. code-block:: shell

	docker-machine start default
	eval $(docker-machine env default)

Execute
------------
.. code-block:: shell

	./nextflow run nmdp-bioinformatics/flow-blast-hml -with-docker \
	nmdpbioinformatics/docker-blast-hml \
	--hml tutorial/ex00_ngsp_expected.xml \
	--outdir tutorial/output --report 0
