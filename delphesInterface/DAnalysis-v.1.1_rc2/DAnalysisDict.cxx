// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME DAnalysisDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "classes/DelphesClasses.h"
#include "vector"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *vectorlEWeightgR_Dictionary();
   static void vectorlEWeightgR_TClassManip(TClass*);
   static void *new_vectorlEWeightgR(void *p = 0);
   static void *newArray_vectorlEWeightgR(Long_t size, void *p);
   static void delete_vectorlEWeightgR(void *p);
   static void deleteArray_vectorlEWeightgR(void *p);
   static void destruct_vectorlEWeightgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<Weight>*)
   {
      vector<Weight> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<Weight>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<Weight>", -2, "vector", 214,
                  typeid(vector<Weight>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEWeightgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<Weight>) );
      instance.SetNew(&new_vectorlEWeightgR);
      instance.SetNewArray(&newArray_vectorlEWeightgR);
      instance.SetDelete(&delete_vectorlEWeightgR);
      instance.SetDeleteArray(&deleteArray_vectorlEWeightgR);
      instance.SetDestructor(&destruct_vectorlEWeightgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<Weight> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<Weight>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEWeightgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<Weight>*)0x0)->GetClass();
      vectorlEWeightgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEWeightgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEWeightgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Weight> : new vector<Weight>;
   }
   static void *newArray_vectorlEWeightgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Weight>[nElements] : new vector<Weight>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEWeightgR(void *p) {
      delete ((vector<Weight>*)p);
   }
   static void deleteArray_vectorlEWeightgR(void *p) {
      delete [] ((vector<Weight>*)p);
   }
   static void destruct_vectorlEWeightgR(void *p) {
      typedef vector<Weight> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<Weight>

namespace ROOT {
   static TClass *vectorlEVertexgR_Dictionary();
   static void vectorlEVertexgR_TClassManip(TClass*);
   static void *new_vectorlEVertexgR(void *p = 0);
   static void *newArray_vectorlEVertexgR(Long_t size, void *p);
   static void delete_vectorlEVertexgR(void *p);
   static void deleteArray_vectorlEVertexgR(void *p);
   static void destruct_vectorlEVertexgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<Vertex>*)
   {
      vector<Vertex> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<Vertex>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<Vertex>", -2, "vector", 214,
                  typeid(vector<Vertex>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEVertexgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<Vertex>) );
      instance.SetNew(&new_vectorlEVertexgR);
      instance.SetNewArray(&newArray_vectorlEVertexgR);
      instance.SetDelete(&delete_vectorlEVertexgR);
      instance.SetDeleteArray(&deleteArray_vectorlEVertexgR);
      instance.SetDestructor(&destruct_vectorlEVertexgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<Vertex> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<Vertex>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEVertexgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<Vertex>*)0x0)->GetClass();
      vectorlEVertexgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEVertexgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEVertexgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Vertex> : new vector<Vertex>;
   }
   static void *newArray_vectorlEVertexgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Vertex>[nElements] : new vector<Vertex>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEVertexgR(void *p) {
      delete ((vector<Vertex>*)p);
   }
   static void deleteArray_vectorlEVertexgR(void *p) {
      delete [] ((vector<Vertex>*)p);
   }
   static void destruct_vectorlEVertexgR(void *p) {
      typedef vector<Vertex> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<Vertex>

namespace ROOT {
   static TClass *vectorlETrackgR_Dictionary();
   static void vectorlETrackgR_TClassManip(TClass*);
   static void *new_vectorlETrackgR(void *p = 0);
   static void *newArray_vectorlETrackgR(Long_t size, void *p);
   static void delete_vectorlETrackgR(void *p);
   static void deleteArray_vectorlETrackgR(void *p);
   static void destruct_vectorlETrackgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<Track>*)
   {
      vector<Track> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<Track>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<Track>", -2, "vector", 214,
                  typeid(vector<Track>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETrackgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<Track>) );
      instance.SetNew(&new_vectorlETrackgR);
      instance.SetNewArray(&newArray_vectorlETrackgR);
      instance.SetDelete(&delete_vectorlETrackgR);
      instance.SetDeleteArray(&deleteArray_vectorlETrackgR);
      instance.SetDestructor(&destruct_vectorlETrackgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<Track> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<Track>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETrackgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<Track>*)0x0)->GetClass();
      vectorlETrackgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETrackgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETrackgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Track> : new vector<Track>;
   }
   static void *newArray_vectorlETrackgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Track>[nElements] : new vector<Track>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETrackgR(void *p) {
      delete ((vector<Track>*)p);
   }
   static void deleteArray_vectorlETrackgR(void *p) {
      delete [] ((vector<Track>*)p);
   }
   static void destruct_vectorlETrackgR(void *p) {
      typedef vector<Track> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<Track>

namespace ROOT {
   static TClass *vectorlETowergR_Dictionary();
   static void vectorlETowergR_TClassManip(TClass*);
   static void *new_vectorlETowergR(void *p = 0);
   static void *newArray_vectorlETowergR(Long_t size, void *p);
   static void delete_vectorlETowergR(void *p);
   static void deleteArray_vectorlETowergR(void *p);
   static void destruct_vectorlETowergR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<Tower>*)
   {
      vector<Tower> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<Tower>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<Tower>", -2, "vector", 214,
                  typeid(vector<Tower>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETowergR_Dictionary, isa_proxy, 4,
                  sizeof(vector<Tower>) );
      instance.SetNew(&new_vectorlETowergR);
      instance.SetNewArray(&newArray_vectorlETowergR);
      instance.SetDelete(&delete_vectorlETowergR);
      instance.SetDeleteArray(&deleteArray_vectorlETowergR);
      instance.SetDestructor(&destruct_vectorlETowergR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<Tower> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<Tower>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETowergR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<Tower>*)0x0)->GetClass();
      vectorlETowergR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETowergR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETowergR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Tower> : new vector<Tower>;
   }
   static void *newArray_vectorlETowergR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Tower>[nElements] : new vector<Tower>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETowergR(void *p) {
      delete ((vector<Tower>*)p);
   }
   static void deleteArray_vectorlETowergR(void *p) {
      delete [] ((vector<Tower>*)p);
   }
   static void destruct_vectorlETowergR(void *p) {
      typedef vector<Tower> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<Tower>

namespace ROOT {
   static TClass *vectorlEScalarHTgR_Dictionary();
   static void vectorlEScalarHTgR_TClassManip(TClass*);
   static void *new_vectorlEScalarHTgR(void *p = 0);
   static void *newArray_vectorlEScalarHTgR(Long_t size, void *p);
   static void delete_vectorlEScalarHTgR(void *p);
   static void deleteArray_vectorlEScalarHTgR(void *p);
   static void destruct_vectorlEScalarHTgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<ScalarHT>*)
   {
      vector<ScalarHT> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<ScalarHT>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<ScalarHT>", -2, "vector", 214,
                  typeid(vector<ScalarHT>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEScalarHTgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<ScalarHT>) );
      instance.SetNew(&new_vectorlEScalarHTgR);
      instance.SetNewArray(&newArray_vectorlEScalarHTgR);
      instance.SetDelete(&delete_vectorlEScalarHTgR);
      instance.SetDeleteArray(&deleteArray_vectorlEScalarHTgR);
      instance.SetDestructor(&destruct_vectorlEScalarHTgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<ScalarHT> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<ScalarHT>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEScalarHTgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<ScalarHT>*)0x0)->GetClass();
      vectorlEScalarHTgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEScalarHTgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEScalarHTgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<ScalarHT> : new vector<ScalarHT>;
   }
   static void *newArray_vectorlEScalarHTgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<ScalarHT>[nElements] : new vector<ScalarHT>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEScalarHTgR(void *p) {
      delete ((vector<ScalarHT>*)p);
   }
   static void deleteArray_vectorlEScalarHTgR(void *p) {
      delete [] ((vector<ScalarHT>*)p);
   }
   static void destruct_vectorlEScalarHTgR(void *p) {
      typedef vector<ScalarHT> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<ScalarHT>

namespace ROOT {
   static TClass *vectorlERhogR_Dictionary();
   static void vectorlERhogR_TClassManip(TClass*);
   static void *new_vectorlERhogR(void *p = 0);
   static void *newArray_vectorlERhogR(Long_t size, void *p);
   static void delete_vectorlERhogR(void *p);
   static void deleteArray_vectorlERhogR(void *p);
   static void destruct_vectorlERhogR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<Rho>*)
   {
      vector<Rho> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<Rho>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<Rho>", -2, "vector", 214,
                  typeid(vector<Rho>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlERhogR_Dictionary, isa_proxy, 4,
                  sizeof(vector<Rho>) );
      instance.SetNew(&new_vectorlERhogR);
      instance.SetNewArray(&newArray_vectorlERhogR);
      instance.SetDelete(&delete_vectorlERhogR);
      instance.SetDeleteArray(&deleteArray_vectorlERhogR);
      instance.SetDestructor(&destruct_vectorlERhogR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<Rho> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<Rho>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlERhogR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<Rho>*)0x0)->GetClass();
      vectorlERhogR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlERhogR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlERhogR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Rho> : new vector<Rho>;
   }
   static void *newArray_vectorlERhogR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Rho>[nElements] : new vector<Rho>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlERhogR(void *p) {
      delete ((vector<Rho>*)p);
   }
   static void deleteArray_vectorlERhogR(void *p) {
      delete [] ((vector<Rho>*)p);
   }
   static void destruct_vectorlERhogR(void *p) {
      typedef vector<Rho> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<Rho>

namespace ROOT {
   static TClass *vectorlEPhotongR_Dictionary();
   static void vectorlEPhotongR_TClassManip(TClass*);
   static void *new_vectorlEPhotongR(void *p = 0);
   static void *newArray_vectorlEPhotongR(Long_t size, void *p);
   static void delete_vectorlEPhotongR(void *p);
   static void deleteArray_vectorlEPhotongR(void *p);
   static void destruct_vectorlEPhotongR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<Photon>*)
   {
      vector<Photon> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<Photon>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<Photon>", -2, "vector", 214,
                  typeid(vector<Photon>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEPhotongR_Dictionary, isa_proxy, 4,
                  sizeof(vector<Photon>) );
      instance.SetNew(&new_vectorlEPhotongR);
      instance.SetNewArray(&newArray_vectorlEPhotongR);
      instance.SetDelete(&delete_vectorlEPhotongR);
      instance.SetDeleteArray(&deleteArray_vectorlEPhotongR);
      instance.SetDestructor(&destruct_vectorlEPhotongR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<Photon> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<Photon>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEPhotongR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<Photon>*)0x0)->GetClass();
      vectorlEPhotongR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEPhotongR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEPhotongR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Photon> : new vector<Photon>;
   }
   static void *newArray_vectorlEPhotongR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Photon>[nElements] : new vector<Photon>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEPhotongR(void *p) {
      delete ((vector<Photon>*)p);
   }
   static void deleteArray_vectorlEPhotongR(void *p) {
      delete [] ((vector<Photon>*)p);
   }
   static void destruct_vectorlEPhotongR(void *p) {
      typedef vector<Photon> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<Photon>

namespace ROOT {
   static TClass *vectorlEMuongR_Dictionary();
   static void vectorlEMuongR_TClassManip(TClass*);
   static void *new_vectorlEMuongR(void *p = 0);
   static void *newArray_vectorlEMuongR(Long_t size, void *p);
   static void delete_vectorlEMuongR(void *p);
   static void deleteArray_vectorlEMuongR(void *p);
   static void destruct_vectorlEMuongR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<Muon>*)
   {
      vector<Muon> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<Muon>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<Muon>", -2, "vector", 214,
                  typeid(vector<Muon>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEMuongR_Dictionary, isa_proxy, 4,
                  sizeof(vector<Muon>) );
      instance.SetNew(&new_vectorlEMuongR);
      instance.SetNewArray(&newArray_vectorlEMuongR);
      instance.SetDelete(&delete_vectorlEMuongR);
      instance.SetDeleteArray(&deleteArray_vectorlEMuongR);
      instance.SetDestructor(&destruct_vectorlEMuongR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<Muon> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<Muon>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEMuongR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<Muon>*)0x0)->GetClass();
      vectorlEMuongR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEMuongR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEMuongR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Muon> : new vector<Muon>;
   }
   static void *newArray_vectorlEMuongR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Muon>[nElements] : new vector<Muon>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEMuongR(void *p) {
      delete ((vector<Muon>*)p);
   }
   static void deleteArray_vectorlEMuongR(void *p) {
      delete [] ((vector<Muon>*)p);
   }
   static void destruct_vectorlEMuongR(void *p) {
      typedef vector<Muon> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<Muon>

namespace ROOT {
   static TClass *vectorlEMissingETgR_Dictionary();
   static void vectorlEMissingETgR_TClassManip(TClass*);
   static void *new_vectorlEMissingETgR(void *p = 0);
   static void *newArray_vectorlEMissingETgR(Long_t size, void *p);
   static void delete_vectorlEMissingETgR(void *p);
   static void deleteArray_vectorlEMissingETgR(void *p);
   static void destruct_vectorlEMissingETgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<MissingET>*)
   {
      vector<MissingET> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<MissingET>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<MissingET>", -2, "vector", 214,
                  typeid(vector<MissingET>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEMissingETgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<MissingET>) );
      instance.SetNew(&new_vectorlEMissingETgR);
      instance.SetNewArray(&newArray_vectorlEMissingETgR);
      instance.SetDelete(&delete_vectorlEMissingETgR);
      instance.SetDeleteArray(&deleteArray_vectorlEMissingETgR);
      instance.SetDestructor(&destruct_vectorlEMissingETgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<MissingET> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<MissingET>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEMissingETgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<MissingET>*)0x0)->GetClass();
      vectorlEMissingETgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEMissingETgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEMissingETgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<MissingET> : new vector<MissingET>;
   }
   static void *newArray_vectorlEMissingETgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<MissingET>[nElements] : new vector<MissingET>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEMissingETgR(void *p) {
      delete ((vector<MissingET>*)p);
   }
   static void deleteArray_vectorlEMissingETgR(void *p) {
      delete [] ((vector<MissingET>*)p);
   }
   static void destruct_vectorlEMissingETgR(void *p) {
      typedef vector<MissingET> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<MissingET>

namespace ROOT {
   static TClass *vectorlELHEFWeightgR_Dictionary();
   static void vectorlELHEFWeightgR_TClassManip(TClass*);
   static void *new_vectorlELHEFWeightgR(void *p = 0);
   static void *newArray_vectorlELHEFWeightgR(Long_t size, void *p);
   static void delete_vectorlELHEFWeightgR(void *p);
   static void deleteArray_vectorlELHEFWeightgR(void *p);
   static void destruct_vectorlELHEFWeightgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<LHEFWeight>*)
   {
      vector<LHEFWeight> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<LHEFWeight>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<LHEFWeight>", -2, "vector", 214,
                  typeid(vector<LHEFWeight>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlELHEFWeightgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<LHEFWeight>) );
      instance.SetNew(&new_vectorlELHEFWeightgR);
      instance.SetNewArray(&newArray_vectorlELHEFWeightgR);
      instance.SetDelete(&delete_vectorlELHEFWeightgR);
      instance.SetDeleteArray(&deleteArray_vectorlELHEFWeightgR);
      instance.SetDestructor(&destruct_vectorlELHEFWeightgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<LHEFWeight> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<LHEFWeight>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlELHEFWeightgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<LHEFWeight>*)0x0)->GetClass();
      vectorlELHEFWeightgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlELHEFWeightgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlELHEFWeightgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<LHEFWeight> : new vector<LHEFWeight>;
   }
   static void *newArray_vectorlELHEFWeightgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<LHEFWeight>[nElements] : new vector<LHEFWeight>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlELHEFWeightgR(void *p) {
      delete ((vector<LHEFWeight>*)p);
   }
   static void deleteArray_vectorlELHEFWeightgR(void *p) {
      delete [] ((vector<LHEFWeight>*)p);
   }
   static void destruct_vectorlELHEFWeightgR(void *p) {
      typedef vector<LHEFWeight> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<LHEFWeight>

namespace ROOT {
   static TClass *vectorlELHEFEventgR_Dictionary();
   static void vectorlELHEFEventgR_TClassManip(TClass*);
   static void *new_vectorlELHEFEventgR(void *p = 0);
   static void *newArray_vectorlELHEFEventgR(Long_t size, void *p);
   static void delete_vectorlELHEFEventgR(void *p);
   static void deleteArray_vectorlELHEFEventgR(void *p);
   static void destruct_vectorlELHEFEventgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<LHEFEvent>*)
   {
      vector<LHEFEvent> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<LHEFEvent>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<LHEFEvent>", -2, "vector", 214,
                  typeid(vector<LHEFEvent>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlELHEFEventgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<LHEFEvent>) );
      instance.SetNew(&new_vectorlELHEFEventgR);
      instance.SetNewArray(&newArray_vectorlELHEFEventgR);
      instance.SetDelete(&delete_vectorlELHEFEventgR);
      instance.SetDeleteArray(&deleteArray_vectorlELHEFEventgR);
      instance.SetDestructor(&destruct_vectorlELHEFEventgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<LHEFEvent> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<LHEFEvent>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlELHEFEventgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<LHEFEvent>*)0x0)->GetClass();
      vectorlELHEFEventgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlELHEFEventgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlELHEFEventgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<LHEFEvent> : new vector<LHEFEvent>;
   }
   static void *newArray_vectorlELHEFEventgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<LHEFEvent>[nElements] : new vector<LHEFEvent>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlELHEFEventgR(void *p) {
      delete ((vector<LHEFEvent>*)p);
   }
   static void deleteArray_vectorlELHEFEventgR(void *p) {
      delete [] ((vector<LHEFEvent>*)p);
   }
   static void destruct_vectorlELHEFEventgR(void *p) {
      typedef vector<LHEFEvent> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<LHEFEvent>

namespace ROOT {
   static TClass *vectorlELHCOEventgR_Dictionary();
   static void vectorlELHCOEventgR_TClassManip(TClass*);
   static void *new_vectorlELHCOEventgR(void *p = 0);
   static void *newArray_vectorlELHCOEventgR(Long_t size, void *p);
   static void delete_vectorlELHCOEventgR(void *p);
   static void deleteArray_vectorlELHCOEventgR(void *p);
   static void destruct_vectorlELHCOEventgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<LHCOEvent>*)
   {
      vector<LHCOEvent> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<LHCOEvent>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<LHCOEvent>", -2, "vector", 214,
                  typeid(vector<LHCOEvent>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlELHCOEventgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<LHCOEvent>) );
      instance.SetNew(&new_vectorlELHCOEventgR);
      instance.SetNewArray(&newArray_vectorlELHCOEventgR);
      instance.SetDelete(&delete_vectorlELHCOEventgR);
      instance.SetDeleteArray(&deleteArray_vectorlELHCOEventgR);
      instance.SetDestructor(&destruct_vectorlELHCOEventgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<LHCOEvent> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<LHCOEvent>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlELHCOEventgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<LHCOEvent>*)0x0)->GetClass();
      vectorlELHCOEventgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlELHCOEventgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlELHCOEventgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<LHCOEvent> : new vector<LHCOEvent>;
   }
   static void *newArray_vectorlELHCOEventgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<LHCOEvent>[nElements] : new vector<LHCOEvent>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlELHCOEventgR(void *p) {
      delete ((vector<LHCOEvent>*)p);
   }
   static void deleteArray_vectorlELHCOEventgR(void *p) {
      delete [] ((vector<LHCOEvent>*)p);
   }
   static void destruct_vectorlELHCOEventgR(void *p) {
      typedef vector<LHCOEvent> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<LHCOEvent>

namespace ROOT {
   static TClass *vectorlEJetgR_Dictionary();
   static void vectorlEJetgR_TClassManip(TClass*);
   static void *new_vectorlEJetgR(void *p = 0);
   static void *newArray_vectorlEJetgR(Long_t size, void *p);
   static void delete_vectorlEJetgR(void *p);
   static void deleteArray_vectorlEJetgR(void *p);
   static void destruct_vectorlEJetgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<Jet>*)
   {
      vector<Jet> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<Jet>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<Jet>", -2, "vector", 214,
                  typeid(vector<Jet>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEJetgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<Jet>) );
      instance.SetNew(&new_vectorlEJetgR);
      instance.SetNewArray(&newArray_vectorlEJetgR);
      instance.SetDelete(&delete_vectorlEJetgR);
      instance.SetDeleteArray(&deleteArray_vectorlEJetgR);
      instance.SetDestructor(&destruct_vectorlEJetgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<Jet> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<Jet>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEJetgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<Jet>*)0x0)->GetClass();
      vectorlEJetgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEJetgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEJetgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Jet> : new vector<Jet>;
   }
   static void *newArray_vectorlEJetgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Jet>[nElements] : new vector<Jet>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEJetgR(void *p) {
      delete ((vector<Jet>*)p);
   }
   static void deleteArray_vectorlEJetgR(void *p) {
      delete [] ((vector<Jet>*)p);
   }
   static void destruct_vectorlEJetgR(void *p) {
      typedef vector<Jet> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<Jet>

namespace ROOT {
   static TClass *vectorlEHepMCEventgR_Dictionary();
   static void vectorlEHepMCEventgR_TClassManip(TClass*);
   static void *new_vectorlEHepMCEventgR(void *p = 0);
   static void *newArray_vectorlEHepMCEventgR(Long_t size, void *p);
   static void delete_vectorlEHepMCEventgR(void *p);
   static void deleteArray_vectorlEHepMCEventgR(void *p);
   static void destruct_vectorlEHepMCEventgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<HepMCEvent>*)
   {
      vector<HepMCEvent> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<HepMCEvent>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<HepMCEvent>", -2, "vector", 214,
                  typeid(vector<HepMCEvent>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEHepMCEventgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<HepMCEvent>) );
      instance.SetNew(&new_vectorlEHepMCEventgR);
      instance.SetNewArray(&newArray_vectorlEHepMCEventgR);
      instance.SetDelete(&delete_vectorlEHepMCEventgR);
      instance.SetDeleteArray(&deleteArray_vectorlEHepMCEventgR);
      instance.SetDestructor(&destruct_vectorlEHepMCEventgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<HepMCEvent> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<HepMCEvent>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEHepMCEventgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<HepMCEvent>*)0x0)->GetClass();
      vectorlEHepMCEventgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEHepMCEventgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEHepMCEventgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<HepMCEvent> : new vector<HepMCEvent>;
   }
   static void *newArray_vectorlEHepMCEventgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<HepMCEvent>[nElements] : new vector<HepMCEvent>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEHepMCEventgR(void *p) {
      delete ((vector<HepMCEvent>*)p);
   }
   static void deleteArray_vectorlEHepMCEventgR(void *p) {
      delete [] ((vector<HepMCEvent>*)p);
   }
   static void destruct_vectorlEHepMCEventgR(void *p) {
      typedef vector<HepMCEvent> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<HepMCEvent>

namespace ROOT {
   static TClass *vectorlEHectorHitgR_Dictionary();
   static void vectorlEHectorHitgR_TClassManip(TClass*);
   static void *new_vectorlEHectorHitgR(void *p = 0);
   static void *newArray_vectorlEHectorHitgR(Long_t size, void *p);
   static void delete_vectorlEHectorHitgR(void *p);
   static void deleteArray_vectorlEHectorHitgR(void *p);
   static void destruct_vectorlEHectorHitgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<HectorHit>*)
   {
      vector<HectorHit> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<HectorHit>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<HectorHit>", -2, "vector", 214,
                  typeid(vector<HectorHit>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEHectorHitgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<HectorHit>) );
      instance.SetNew(&new_vectorlEHectorHitgR);
      instance.SetNewArray(&newArray_vectorlEHectorHitgR);
      instance.SetDelete(&delete_vectorlEHectorHitgR);
      instance.SetDeleteArray(&deleteArray_vectorlEHectorHitgR);
      instance.SetDestructor(&destruct_vectorlEHectorHitgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<HectorHit> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<HectorHit>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEHectorHitgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<HectorHit>*)0x0)->GetClass();
      vectorlEHectorHitgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEHectorHitgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEHectorHitgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<HectorHit> : new vector<HectorHit>;
   }
   static void *newArray_vectorlEHectorHitgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<HectorHit>[nElements] : new vector<HectorHit>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEHectorHitgR(void *p) {
      delete ((vector<HectorHit>*)p);
   }
   static void deleteArray_vectorlEHectorHitgR(void *p) {
      delete [] ((vector<HectorHit>*)p);
   }
   static void destruct_vectorlEHectorHitgR(void *p) {
      typedef vector<HectorHit> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<HectorHit>

namespace ROOT {
   static TClass *vectorlEGenParticlegR_Dictionary();
   static void vectorlEGenParticlegR_TClassManip(TClass*);
   static void *new_vectorlEGenParticlegR(void *p = 0);
   static void *newArray_vectorlEGenParticlegR(Long_t size, void *p);
   static void delete_vectorlEGenParticlegR(void *p);
   static void deleteArray_vectorlEGenParticlegR(void *p);
   static void destruct_vectorlEGenParticlegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<GenParticle>*)
   {
      vector<GenParticle> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<GenParticle>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<GenParticle>", -2, "vector", 214,
                  typeid(vector<GenParticle>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEGenParticlegR_Dictionary, isa_proxy, 4,
                  sizeof(vector<GenParticle>) );
      instance.SetNew(&new_vectorlEGenParticlegR);
      instance.SetNewArray(&newArray_vectorlEGenParticlegR);
      instance.SetDelete(&delete_vectorlEGenParticlegR);
      instance.SetDeleteArray(&deleteArray_vectorlEGenParticlegR);
      instance.SetDestructor(&destruct_vectorlEGenParticlegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<GenParticle> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<GenParticle>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEGenParticlegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<GenParticle>*)0x0)->GetClass();
      vectorlEGenParticlegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEGenParticlegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEGenParticlegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<GenParticle> : new vector<GenParticle>;
   }
   static void *newArray_vectorlEGenParticlegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<GenParticle>[nElements] : new vector<GenParticle>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEGenParticlegR(void *p) {
      delete ((vector<GenParticle>*)p);
   }
   static void deleteArray_vectorlEGenParticlegR(void *p) {
      delete [] ((vector<GenParticle>*)p);
   }
   static void destruct_vectorlEGenParticlegR(void *p) {
      typedef vector<GenParticle> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<GenParticle>

namespace ROOT {
   static TClass *vectorlEEventgR_Dictionary();
   static void vectorlEEventgR_TClassManip(TClass*);
   static void *new_vectorlEEventgR(void *p = 0);
   static void *newArray_vectorlEEventgR(Long_t size, void *p);
   static void delete_vectorlEEventgR(void *p);
   static void deleteArray_vectorlEEventgR(void *p);
   static void destruct_vectorlEEventgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<Event>*)
   {
      vector<Event> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<Event>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<Event>", -2, "vector", 214,
                  typeid(vector<Event>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEEventgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<Event>) );
      instance.SetNew(&new_vectorlEEventgR);
      instance.SetNewArray(&newArray_vectorlEEventgR);
      instance.SetDelete(&delete_vectorlEEventgR);
      instance.SetDeleteArray(&deleteArray_vectorlEEventgR);
      instance.SetDestructor(&destruct_vectorlEEventgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<Event> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<Event>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEEventgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<Event>*)0x0)->GetClass();
      vectorlEEventgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEEventgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEEventgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Event> : new vector<Event>;
   }
   static void *newArray_vectorlEEventgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Event>[nElements] : new vector<Event>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEEventgR(void *p) {
      delete ((vector<Event>*)p);
   }
   static void deleteArray_vectorlEEventgR(void *p) {
      delete [] ((vector<Event>*)p);
   }
   static void destruct_vectorlEEventgR(void *p) {
      typedef vector<Event> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<Event>

namespace ROOT {
   static TClass *vectorlEElectrongR_Dictionary();
   static void vectorlEElectrongR_TClassManip(TClass*);
   static void *new_vectorlEElectrongR(void *p = 0);
   static void *newArray_vectorlEElectrongR(Long_t size, void *p);
   static void delete_vectorlEElectrongR(void *p);
   static void deleteArray_vectorlEElectrongR(void *p);
   static void destruct_vectorlEElectrongR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<Electron>*)
   {
      vector<Electron> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<Electron>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<Electron>", -2, "vector", 214,
                  typeid(vector<Electron>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEElectrongR_Dictionary, isa_proxy, 4,
                  sizeof(vector<Electron>) );
      instance.SetNew(&new_vectorlEElectrongR);
      instance.SetNewArray(&newArray_vectorlEElectrongR);
      instance.SetDelete(&delete_vectorlEElectrongR);
      instance.SetDeleteArray(&deleteArray_vectorlEElectrongR);
      instance.SetDestructor(&destruct_vectorlEElectrongR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<Electron> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<Electron>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEElectrongR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<Electron>*)0x0)->GetClass();
      vectorlEElectrongR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEElectrongR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEElectrongR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Electron> : new vector<Electron>;
   }
   static void *newArray_vectorlEElectrongR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Electron>[nElements] : new vector<Electron>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEElectrongR(void *p) {
      delete ((vector<Electron>*)p);
   }
   static void deleteArray_vectorlEElectrongR(void *p) {
      delete [] ((vector<Electron>*)p);
   }
   static void destruct_vectorlEElectrongR(void *p) {
      typedef vector<Electron> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<Electron>

namespace ROOT {
   static TClass *vectorlECandidategR_Dictionary();
   static void vectorlECandidategR_TClassManip(TClass*);
   static void *new_vectorlECandidategR(void *p = 0);
   static void *newArray_vectorlECandidategR(Long_t size, void *p);
   static void delete_vectorlECandidategR(void *p);
   static void deleteArray_vectorlECandidategR(void *p);
   static void destruct_vectorlECandidategR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<Candidate>*)
   {
      vector<Candidate> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<Candidate>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<Candidate>", -2, "vector", 214,
                  typeid(vector<Candidate>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlECandidategR_Dictionary, isa_proxy, 4,
                  sizeof(vector<Candidate>) );
      instance.SetNew(&new_vectorlECandidategR);
      instance.SetNewArray(&newArray_vectorlECandidategR);
      instance.SetDelete(&delete_vectorlECandidategR);
      instance.SetDeleteArray(&deleteArray_vectorlECandidategR);
      instance.SetDestructor(&destruct_vectorlECandidategR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<Candidate> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<Candidate>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlECandidategR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<Candidate>*)0x0)->GetClass();
      vectorlECandidategR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlECandidategR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlECandidategR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Candidate> : new vector<Candidate>;
   }
   static void *newArray_vectorlECandidategR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Candidate>[nElements] : new vector<Candidate>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlECandidategR(void *p) {
      delete ((vector<Candidate>*)p);
   }
   static void deleteArray_vectorlECandidategR(void *p) {
      delete [] ((vector<Candidate>*)p);
   }
   static void destruct_vectorlECandidategR(void *p) {
      typedef vector<Candidate> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<Candidate>

namespace {
  void TriggerDictionaryInitialization_DAnalysisDict_Impl() {
    static const char* headers[] = {
"classes/DelphesClasses.h",
"vector",
0
    };
    static const char* includePaths[] = {
"/cvmfs/cms.cern.ch/slc6_amd64_gcc630/cms/cmssw/CMSSW_9_3_2/external/slc6_amd64_gcc630/bin/../../../../../../../slc6_amd64_gcc630/lcg/root/6.10.04-ghjeda/include",
"/afs/cern.ch/work/j/jkiesele/Phase2/CMSSW_9_3_2/src/PhaseTwoAnalysis/delphesInterface/delphes",
"/cvmfs/cms.cern.ch/slc6_amd64_gcc630/lcg/root/6.10.04-ghjeda/include",
"/afs/cern.ch/work/j/jkiesele/Phase2/CMSSW_9_3_2/src/PhaseTwoAnalysis/delphesInterface/DAnalysis-v.1.1_rc2/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "DAnalysisDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$classes/DelphesClasses.h")))  Candidate;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
class __attribute__((annotate("$clingAutoload$classes/DelphesClasses.h")))  HectorHit;
class __attribute__((annotate("$clingAutoload$classes/DelphesClasses.h")))  Tower;
class __attribute__((annotate("$clingAutoload$classes/DelphesClasses.h")))  Track;
class __attribute__((annotate("$clingAutoload$classes/DelphesClasses.h")))  Jet;
class __attribute__((annotate("$clingAutoload$classes/DelphesClasses.h")))  Muon;
class __attribute__((annotate("$clingAutoload$classes/DelphesClasses.h")))  Electron;
class __attribute__((annotate("$clingAutoload$classes/DelphesClasses.h")))  Photon;
class __attribute__((annotate("$clingAutoload$classes/DelphesClasses.h")))  Weight;
class __attribute__((annotate("$clingAutoload$classes/DelphesClasses.h")))  Rho;
class __attribute__((annotate("$clingAutoload$classes/DelphesClasses.h")))  ScalarHT;
class __attribute__((annotate("$clingAutoload$classes/DelphesClasses.h")))  MissingET;
class __attribute__((annotate("$clingAutoload$classes/DelphesClasses.h")))  Vertex;
class __attribute__((annotate("$clingAutoload$classes/DelphesClasses.h")))  GenParticle;
class __attribute__((annotate("$clingAutoload$classes/DelphesClasses.h")))  HepMCEvent;
class __attribute__((annotate("$clingAutoload$classes/DelphesClasses.h")))  LHEFWeight;
class __attribute__((annotate("$clingAutoload$classes/DelphesClasses.h")))  LHEFEvent;
class __attribute__((annotate("$clingAutoload$classes/DelphesClasses.h")))  LHCOEvent;
class __attribute__((annotate("$clingAutoload$classes/DelphesClasses.h")))  Event;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "DAnalysisDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "classes/DelphesClasses.h"
#include "vector"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("DAnalysisDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_DAnalysisDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_DAnalysisDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_DAnalysisDict() {
  TriggerDictionaryInitialization_DAnalysisDict_Impl();
}
