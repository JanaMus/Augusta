Installation
------------

Augusta is an open-source **Python 3 (up to 3.8)** package.
The lines bellow provides through the installation.

Dependencies: docker

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


For further information see `docs.docker.com/engine/install/linux-postinstall <https://docs.docker.com/engine/install/linux-postinstall/>`_.
