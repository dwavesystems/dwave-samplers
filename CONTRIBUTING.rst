Contributing
============

Ocean's `contributing guide <https://docs.dwavequantum.com/en/latest/ocean/contribute.html>`_
has guidelines for contributing to Ocean packages.

Release Notes
-------------

**dwave-samplers** makes use of `reno <https://docs.openstack.org/reno/>`_ to
manage its release notes.

When making a contribution to **dwave-samplers** that will affect users, create
a new release note file by running

.. code-block:: bash

    reno new your-short-descriptor-here

You can then edit the file created under ``releasenotes/notes/``.
Remove any sections not relevant to your changes.
Commit the file along with your changes.
