.. _cpp_docs:

C++ API Documentation
---------------------

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

.. doxygenstruct:: grgl::NodeData
    :members:


Serialization and Conversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: grgl::writeGrg

.. doxygenfunction:: grgl::readMutableGrg

.. doxygenfunction:: grgl::readImmutableGrg

.. doxygenfunction:: grgl::convertTreeSeqToGRG

Visitors
~~~~~~~~

.. doxygenenum:: grgl::TraversalDirection

.. doxygenenum:: grgl::DfsPass

.. doxygenclass:: grgl::GRGVisitor
    :members:

