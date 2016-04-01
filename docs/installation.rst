Installation
=============

Git
------------------------

Install nextflow by running the following command.

.. code-block:: shell

	apt-get install git

For further documentation on git `click here <https://git-scm.com/>`_.

Nextflow
------------------------

Install nextflow by running the following command.

.. code-block:: shell

	curl -fsSL get.nextflow.io | bash

For further documentation on nextflow `click here <http://www.nextflow.io/>`_.

Docker
------------------------

`Click here <https://docs.docker.com/>`_ for instructions on how to install docker. After you have docker installed run the following commands to start the docker daemon on your machine.

.. code-block:: shell

	docker-machine start default
	eval $(docker-machine env default)

.. note:: If you're using a VPN the docker-machine will not work.


