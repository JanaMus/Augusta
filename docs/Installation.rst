Installation
------------

We highly reccoment installing and using Augusta in a virtual environment. 

Virtual environment can be created and activated using  `conda <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_:

.. code-block:: 

   > conda create -n venv_Augusta python=3.7 anaconda
   > conda activate venv_Augusta
   

Dependencies
=====================

Docker
^^^^^^^^
Docker installation can be checked via:

.. code-block:: 

   > docker info
   
See  `Docker <https://docs.docker.com/get-docker/>`_ for more information.


Python
^^^^^^^^^
Python is already installed in the virtual environment venv_Augusta. 
In case of not using venv_Augusta, Augusta package needs **Python version 3.7 or 3.8**. 

Python version can be checked via:

.. code-block:: 

   > python3 --version
   
See `Python <https://www.python.org/>`_ for more information.


Install Augusta
==================

from PyPi / pip
^^^^^^^^^^^^^^^^

.. code-block:: python

   > pip install Augusta


from GitHub
^^^^^^^^^^^

.. code-block:: python

   > git clone https://github.com/JanaMus/Augusta.git
   > cd Augusta
   > python setup.py install


*Linux non-root user:*
Augusta uses MEME Suite **Docker** Image to search for motifs.
Therefore, access to Docker is needed to be set via terminal:

1. create the docker group

.. code-block:: python

   > sudo groupadd docker


2. add your user to the docker group

.. code-block:: python

   > sudo usermod -aG docker $USER


See `Docker documentation <https://docs.docker.com/engine/install/linux-postinstall/>`_ for more information.
