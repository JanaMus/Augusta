Installation
------------

Augusta is an open-source, OS independent Python 3 package.
The lines bellow provides through the installation.

from PyPi / pip
^^^^^^^^^^^^^^^^
.. code-block::

   $ pip install Augusta


from GitHub
^^^^^^^^^^^
.. code-block::

   $ git clone https://github.com/JanaMus/Augusta.git
   $ cd Augusta
   $ python setup.py install


*Linux non-root user:*
Augusta uses MEME Suite **Docker** Image to search for motifs.
Therefore, access to Docker is needed to be set via terminal:

1. create the docker group
.. code-block::

   $ sudo groupadd docker


2. add your user to the docker group
.. code-block::

   $ sudo usermod -aG docker $USER


For further information see `docs.docker.com/engine/install/linux-postinstall <https://docs.docker.com/engine/install/linux-postinstall/>`_