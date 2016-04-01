Running the Tutorial
================================

Step 1: Clone github repository
-----------------------------
.. code-block:: shell

	git clone https://github.com/nmdp-bioinformatics/flow-blast-hml
	cd flow-blast-hml

Step 2: Install nextflow
-----------------------------
.. code-block:: shell

	curl -fsSL get.nextflow.io | bash

Step 3: Start docker machine
-----------------------------
.. code-block:: shell

	docker-machine start default
	eval $(docker-machine env default)

Step 4: Execute
-----------------------------
.. code-block:: shell

	./nextflow run nmdp-bioinformatics/flow-blast-hml -with-docker \
	nmdpbioinformatics/docker-blast-hml \
	--hml tutorial/ex00_ngsp_expected.xml \
	--outdir tutorial/output --report 0
