.. _cpp_docs:

C++ API
-------

Exception Types
~~~~~~~~~~~~~~~

.. doxygenclass:: grgl::ApiMisuseFailure
    :members:

.. doxygenclass:: grgl::BadInputFileFailure
    :members:

.. doxygenclass:: grgl::SerializationFailure
    :members:

Main GRG Classes
~~~~~~~~~~~~~~~~

.. doxygenclass:: grgl::GRG
    :members:

.. doxygenclass:: grgl::MutableGRG
    :members:

.. doxygenclass:: grgl::Mutation
    :members:

.. doxygenclass:: grgl::GRG::NodeAndMut
.. doxygenclass:: grgl::GRG::NodeMutMiss
.. doxygenclass:: grgl::GRG::MutAndNode
.. doxygenclass:: grgl::GRG::MutNodeMiss


Serialization and Conversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: grgl::writeGrg

.. doxygenfunction:: grgl::readMutableGrg

.. doxygenfunction:: grgl::readImmutableGrg

.. doxygenfunction:: grgl::convertTreeSeqToGRG

.. doxygenfunction:: grgl::simplifyAndSerialize

Visitors
~~~~~~~~

.. doxygenenum:: grgl::TraversalDirection

.. doxygenenum:: grgl::DfsPass

.. doxygenclass:: grgl::GRGVisitor
    :members:

