Installation
------------

Install dependencies
=====================

Python
^^^^^^^^^
Augusta is an open-source **Python 3 (up to 3.8)** package. 
Python comes preinstalled on most Linux and with OSX distributions. 
On Windows, Python must be installed manually: see `Python´s Documentation <https://docs.python.org/3/using/windows.html>`_.

Python version can be checked via (Augusta works with Python 3, up to 3.8):

.. code-block:: 

   > python3
   
See `Python <https://www.python.org/>`_ for more information.

Docker
^^^^^^^^
Docker installation can be checked via:

.. code-block:: 

   > docker info
   
See  `Docker <https://docs.docker.com/get-docker/>`_ for more information.


Install Augusta
==================

from PyPi / pip
^^^^^^^^^^^^^^^^

.. code-block:: python

   > pip3 install Augusta


from GitHub
^^^^^^^^^^^

.. code-block:: python

   > git clone https://github.com/JanaMus/Augusta.git
   > cd Augusta
   > python3 setup.py install


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
