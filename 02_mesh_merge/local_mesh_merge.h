#ifndef POISSON_MULTIMESH_LOCAL_MESH_MERGE_H
#define POISSON_MULTIMESH_LOCAL_MESH_MERGE_H

#include <cstddef>
#include <stdexcept>

#include <AFEPack/HGeometry.h>

namespace local_multimesh {

// The common mesh is the union of the refinement trees of two IrregularMesh
// objects. Both inputs must have the same root-element ordering and must use
// AFEPack's standard child ordering.
template <int DIM, int DOW = DIM>
class CommonIrregularMesh : public IrregularMesh<DIM, DOW> {
 public:
  using Base = IrregularMesh<DIM, DOW>;
  using Element = HElement<DIM, DOW>;

  explicit CommonIrregularMesh(const Base& first) : Base(first) {}

  void merge(const Base& second) {
    std::size_t second_root_count = 0;
    for (typename Base::ConstRootIterator it = second.beginRootElement();
         it != second.endRootElement(); ++it) {
      ++second_root_count;
    }
    if (this->rootElement().size() != second_root_count) {
      throw std::invalid_argument(
          "Cannot merge meshes with different root-element counts.");
    }

    typename Base::RootIterator destination = this->beginRootElement();
    typename Base::ConstRootIterator source = second.beginRootElement();
    for (; destination != this->endRootElement(); ++destination, ++source) {
      merge_subtree(*destination, *source);
    }
  }

 private:
  void merge_subtree(Element& destination, const Element& source) {
    if (!source.isRefined()) {
      return;
    }

    if (!destination.isRefined()) {
      // refineElement(), rather than HElement::refine(), also sets AFEPack's
      // active-element values on the new children.
      this->refineElement(destination);
    }

    if (destination.child.size() != source.child.size()) {
      throw std::runtime_error(
          "The two meshes use incompatible refinement child layouts.");
    }

    for (std::size_t i = 0; i < source.child.size(); ++i) {
      merge_subtree(*destination.child[i], *source.child[i]);
    }
  }
};

}  // namespace local_multimesh

#endif  // POISSON_MULTIMESH_LOCAL_MESH_MERGE_H
