#include<memory>
#include<vector>

using namespace std;
//container for managing pointers to a class.
//class must implement the method T.deepcopy(), which returns
//a unique_ptr to a copy of the object. Useful for class
//polymorphism - I can make a ptr_vec<BaseClass> and fill it
//with unique_ptrs to DerivedClass objects.

//Richard Taylor, Geoscience Australia 2020.

template<typename T>
class ptr_vec {
private:
  vector<unique_ptr<T>> v;
public:
  ptr_vec() = default;
  ptr_vec(const ptr_vec& other) {
    for (size_t i = 0; i < other.v.size(); i++) {
      v.push_back(other.v[i]->deepcopy());
    }
  }
  ptr_vec& operator=(const ptr_vec& other) {
    //we own the pointers so the objects they point to
    //are destroyed when we do this
    v.clear();
    //make new unique_ptrs
    for (size_t i = 0; i < other.v.size(); i++) {
      v.push_back(other.v[i]->deepcopy());
    }
    return *this;
  }

  ptr_vec(ptr_vec&& other) noexcept {
    v = std::move(other.v);
  }

  ptr_vec& operator=(ptr_vec&& other) noexcept {
    if (this != &other){
      v.clear();
      v = std::move(other.v);
    }
    return *this;
  }

  ~ptr_vec() = default;


  T* operator[](const size_t& idx) const {
    return v[idx].get();
  }

  void push_back(unique_ptr<T> ptr) {
    v.push_back(std::move(ptr));
  }

  void clear() {
    v.clear();
  }

  size_t size() const {
    return v.size();
  }

};
