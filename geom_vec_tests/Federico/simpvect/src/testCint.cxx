//
// File generated by rootcint at Thu Jan 31 19:33:45 2013

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME testCint
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "testCint.h"

#include "TClass.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"

// START OF SHADOWS

namespace ROOT {
   namespace Shadow {
   } // of namespace Shadow
} // of namespace ROOT
// END OF SHADOWS

namespace ROOT {
   void TGeoBBox_v_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_TGeoBBox_v(void *p = 0);
   static void *newArray_TGeoBBox_v(Long_t size, void *p);
   static void delete_TGeoBBox_v(void *p);
   static void deleteArray_TGeoBBox_v(void *p);
   static void destruct_TGeoBBox_v(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TGeoBBox_v*)
   {
      ::TGeoBBox_v *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TGeoBBox_v >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TGeoBBox_v", ::TGeoBBox_v::Class_Version(), "./TGeoBBOx_v.h", 32,
                  typeid(::TGeoBBox_v), DefineBehavior(ptr, ptr),
                  &::TGeoBBox_v::Dictionary, isa_proxy, 4,
                  sizeof(::TGeoBBox_v) );
      instance.SetNew(&new_TGeoBBox_v);
      instance.SetNewArray(&newArray_TGeoBBox_v);
      instance.SetDelete(&delete_TGeoBBox_v);
      instance.SetDeleteArray(&deleteArray_TGeoBBox_v);
      instance.SetDestructor(&destruct_TGeoBBox_v);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TGeoBBox_v*)
   {
      return GenerateInitInstanceLocal((::TGeoBBox_v*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TGeoBBox_v*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
TClass *TGeoBBox_v::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *TGeoBBox_v::Class_Name()
{
   return "TGeoBBox_v";
}

//______________________________________________________________________________
const char *TGeoBBox_v::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TGeoBBox_v*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TGeoBBox_v::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TGeoBBox_v*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void TGeoBBox_v::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TGeoBBox_v*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *TGeoBBox_v::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TGeoBBox_v*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
void TGeoBBox_v::Streamer(TBuffer &R__b)
{
   // Stream an object of class TGeoBBox_v.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TGeoBBox_v::Class(),this);
   } else {
      R__b.WriteClassBuffer(TGeoBBox_v::Class(),this);
   }
}

//______________________________________________________________________________
void TGeoBBox_v::ShowMembers(TMemberInspector &R__insp)
{
      // Inspect the data members of an object of class TGeoBBox_v.
      TClass *R__cl = ::TGeoBBox_v::IsA();
      if (R__cl || R__insp.IsA()) { }
      TGeoBBox::ShowMembers(R__insp);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TGeoBBox_v(void *p) {
      return  p ? new(p) ::TGeoBBox_v : new ::TGeoBBox_v;
   }
   static void *newArray_TGeoBBox_v(Long_t nElements, void *p) {
      return p ? new(p) ::TGeoBBox_v[nElements] : new ::TGeoBBox_v[nElements];
   }
   // Wrapper around operator delete
   static void delete_TGeoBBox_v(void *p) {
      delete ((::TGeoBBox_v*)p);
   }
   static void deleteArray_TGeoBBox_v(void *p) {
      delete [] ((::TGeoBBox_v*)p);
   }
   static void destruct_TGeoBBox_v(void *p) {
      typedef ::TGeoBBox_v current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TGeoBBox_v

/********************************************************
* testCint.cxx
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************/

#ifdef G__MEMTEST
#undef malloc
#undef free
#endif

#if defined(__GNUC__) && __GNUC__ >= 4 && ((__GNUC_MINOR__ == 2 && __GNUC_PATCHLEVEL__ >= 1) || (__GNUC_MINOR__ >= 3))
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

extern "C" void G__cpp_reset_tagtabletestCint();

extern "C" void G__set_cpp_environmenttestCint() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("TGeoBBOx_v.h");
  G__cpp_reset_tagtabletestCint();
}
#include <new>
extern "C" int G__cpp_dllrevtestCint() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* TGeoBBox_v */
static int G__testCint_191_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   TGeoBBox_v* p = NULL;
   char* gvp = (char*) G__getgvp();
   int n = G__getaryconstruct();
   if (n) {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new TGeoBBox_v[n];
     } else {
       p = new((void*) gvp) TGeoBBox_v[n];
     }
   } else {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new TGeoBBox_v;
     } else {
       p = new((void*) gvp) TGeoBBox_v;
     }
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__testCintLN_TGeoBBox_v));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__testCint_191_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   TGeoBBox_v* p = NULL;
   char* gvp = (char*) G__getgvp();
   switch (libp->paran) {
   case 4:
     //m: 4
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new TGeoBBox_v(
(Double_t) G__double(libp->para[0]), (Double_t) G__double(libp->para[1])
, (Double_t) G__double(libp->para[2]), (Double_t*) G__int(libp->para[3]));
     } else {
       p = new((void*) gvp) TGeoBBox_v(
(Double_t) G__double(libp->para[0]), (Double_t) G__double(libp->para[1])
, (Double_t) G__double(libp->para[2]), (Double_t*) G__int(libp->para[3]));
     }
     break;
   case 3:
     //m: 3
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new TGeoBBox_v(
(Double_t) G__double(libp->para[0]), (Double_t) G__double(libp->para[1])
, (Double_t) G__double(libp->para[2]));
     } else {
       p = new((void*) gvp) TGeoBBox_v(
(Double_t) G__double(libp->para[0]), (Double_t) G__double(libp->para[1])
, (Double_t) G__double(libp->para[2]));
     }
     break;
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__testCintLN_TGeoBBox_v));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__testCint_191_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   TGeoBBox_v* p = NULL;
   char* gvp = (char*) G__getgvp();
   switch (libp->paran) {
   case 5:
     //m: 5
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new TGeoBBox_v(
(const char*) G__int(libp->para[0]), (Double_t) G__double(libp->para[1])
, (Double_t) G__double(libp->para[2]), (Double_t) G__double(libp->para[3])
, (Double_t*) G__int(libp->para[4]));
     } else {
       p = new((void*) gvp) TGeoBBox_v(
(const char*) G__int(libp->para[0]), (Double_t) G__double(libp->para[1])
, (Double_t) G__double(libp->para[2]), (Double_t) G__double(libp->para[3])
, (Double_t*) G__int(libp->para[4]));
     }
     break;
   case 4:
     //m: 4
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new TGeoBBox_v(
(const char*) G__int(libp->para[0]), (Double_t) G__double(libp->para[1])
, (Double_t) G__double(libp->para[2]), (Double_t) G__double(libp->para[3]));
     } else {
       p = new((void*) gvp) TGeoBBox_v(
(const char*) G__int(libp->para[0]), (Double_t) G__double(libp->para[1])
, (Double_t) G__double(libp->para[2]), (Double_t) G__double(libp->para[3]));
     }
     break;
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__testCintLN_TGeoBBox_v));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__testCint_191_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   TGeoBBox_v* p = NULL;
   char* gvp = (char*) G__getgvp();
   //m: 1
   if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
     p = new TGeoBBox_v((Double_t*) G__int(libp->para[0]));
   } else {
     p = new((void*) gvp) TGeoBBox_v((Double_t*) G__int(libp->para[0]));
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__testCintLN_TGeoBBox_v));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__testCint_191_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 103, (long) TGeoBBox_v::AreOverlapping((TGeoBBox_v*) G__int(libp->para[0]), (TGeoMatrix*) G__int(libp->para[1])
, (TGeoBBox_v*) G__int(libp->para[2]), (TGeoMatrix*) G__int(libp->para[3])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__testCint_191_0_8(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((const TGeoBBox_v*) G__getstructoffset())->Contains_v((Double_t*) G__int(libp->para[0]), (Bool_t*) G__int(libp->para[1])
, (Int_t) G__int(libp->para[2]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__testCint_191_0_9(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 103, (long) TGeoBBox_v::Contains((Double_t*) G__int(libp->para[0]), (Double_t) G__double(libp->para[1])
, (Double_t) G__double(libp->para[2]), (Double_t) G__double(libp->para[3])
, (Double_t*) G__int(libp->para[4])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__testCint_191_0_13(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   switch (libp->paran) {
   case 7:
      G__letdouble(result7, 100, (double) TGeoBBox_v::DistFromInside(
(Double_t*) G__int(libp->para[0]), (Double_t*) G__int(libp->para[1])
, (Double_t) G__double(libp->para[2]), (Double_t) G__double(libp->para[3])
, (Double_t) G__double(libp->para[4]), (Double_t*) G__int(libp->para[5])
, (Double_t) G__double(libp->para[6])));
      break;
   case 6:
      G__letdouble(result7, 100, (double) TGeoBBox_v::DistFromInside((Double_t*) G__int(libp->para[0]), (Double_t*) G__int(libp->para[1])
, (Double_t) G__double(libp->para[2]), (Double_t) G__double(libp->para[3])
, (Double_t) G__double(libp->para[4]), (Double_t*) G__int(libp->para[5])));
      break;
   }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__testCint_191_0_15(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   switch (libp->paran) {
   case 7:
      G__letdouble(result7, 100, (double) TGeoBBox_v::DistFromOutside(
(Double_t*) G__int(libp->para[0]), (Double_t*) G__int(libp->para[1])
, (Double_t) G__double(libp->para[2]), (Double_t) G__double(libp->para[3])
, (Double_t) G__double(libp->para[4]), (Double_t*) G__int(libp->para[5])
, (Double_t) G__double(libp->para[6])));
      break;
   case 6:
      G__letdouble(result7, 100, (double) TGeoBBox_v::DistFromOutside((Double_t*) G__int(libp->para[0]), (Double_t*) G__int(libp->para[1])
, (Double_t) G__double(libp->para[2]), (Double_t) G__double(libp->para[3])
, (Double_t) G__double(libp->para[4]), (Double_t*) G__int(libp->para[5])));
      break;
   }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__testCint_191_0_19(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) TGeoBBox_v::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__testCint_191_0_20(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) TGeoBBox_v::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__testCint_191_0_21(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) TGeoBBox_v::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__testCint_191_0_22(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      TGeoBBox_v::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__testCint_191_0_26(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((TGeoBBox_v*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__testCint_191_0_27(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) TGeoBBox_v::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__testCint_191_0_28(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) TGeoBBox_v::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__testCint_191_0_29(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) TGeoBBox_v::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__testCint_191_0_30(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) TGeoBBox_v::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__testCint_191_0_31(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   TGeoBBox_v* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new TGeoBBox_v(*(TGeoBBox_v*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__testCintLN_TGeoBBox_v));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef TGeoBBox_v G__TTGeoBBox_v;
static int G__testCint_191_0_32(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   char* gvp = (char*) G__getgvp();
   long soff = G__getstructoffset();
   int n = G__getaryconstruct();
   //
   //has_a_delete: 1
   //has_own_delete1arg: 0
   //has_own_delete2arg: 0
   //
   if (!soff) {
     return(1);
   }
   if (n) {
     if (gvp == (char*)G__PVOID) {
       delete[] (TGeoBBox_v*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((TGeoBBox_v*) (soff+(sizeof(TGeoBBox_v)*i)))->~G__TTGeoBBox_v();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (TGeoBBox_v*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((TGeoBBox_v*) (soff))->~G__TTGeoBBox_v();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__testCint_191_0_33(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   TGeoBBox_v* dest = (TGeoBBox_v*) G__getstructoffset();
   *dest = *(TGeoBBox_v*) libp->para[0].ref;
   const TGeoBBox_v& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* TGeoBBox_v */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfunctestCint {
 public:
  G__Sizep2memfunctestCint(): p(&G__Sizep2memfunctestCint::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfunctestCint::*p)();
};

size_t G__get_sizep2memfunctestCint()
{
  G__Sizep2memfunctestCint a;
  G__setsizep2memfunc((int)a.sizep2memfunc());
  return((size_t)a.sizep2memfunc());
}


/*********************************************************
* virtual base class offset calculation interface
*********************************************************/

   /* Setting up class inheritance */

/*********************************************************
* Inheritance information setup/
*********************************************************/
extern "C" void G__cpp_setup_inheritancetestCint() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__testCintLN_TGeoBBox_v))) {
     TGeoBBox_v *G__Lderived;
     G__Lderived=(TGeoBBox_v*)0x1000;
     {
       TGeoBBox *G__Lpbase=(TGeoBBox*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__testCintLN_TGeoBBox_v),G__get_linked_tagnum(&G__testCintLN_TGeoBBox),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
     {
       TGeoShape *G__Lpbase=(TGeoShape*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__testCintLN_TGeoBBox_v),G__get_linked_tagnum(&G__testCintLN_TGeoShape),(long)G__Lpbase-(long)G__Lderived,1,0);
     }
     {
       TNamed *G__Lpbase=(TNamed*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__testCintLN_TGeoBBox_v),G__get_linked_tagnum(&G__testCintLN_TNamed),(long)G__Lpbase-(long)G__Lderived,1,0);
     }
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__testCintLN_TGeoBBox_v),G__get_linked_tagnum(&G__testCintLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,0);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetabletestCint() {

   /* Setting up typedef entry */
   G__search_typename2("Int_t",105,-1,0,-1);
   G__setnewtype(-1,"Signed integer 4 bytes (int)",0);
   G__search_typename2("Double_t",100,-1,0,-1);
   G__setnewtype(-1,"Double 8 bytes",0);
   G__search_typename2("Bool_t",103,-1,0,-1);
   G__setnewtype(-1,"Boolean (0=false, 1=true) (bool)",0);
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__testCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__testCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__testCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__testCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__testCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__testCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__testCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__testCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__testCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__testCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* TGeoBBox_v */
static void G__setup_memvarTGeoBBox_v(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__testCintLN_TGeoBBox_v));
   { TGeoBBox_v *p; p=(TGeoBBox_v*)0x1000; if (p) { }
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__testCintLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvartestCint() {
}
/***********************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
***********************************************************/

/*********************************************************
* Member function information setup for each class
*********************************************************/
static void G__setup_memfuncTGeoBBox_v(void) {
   /* TGeoBBox_v */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__testCintLN_TGeoBBox_v));
   G__memfunc_setup("TGeoBBox_v",943,G__testCint_191_0_1, 105, G__get_linked_tagnum(&G__testCintLN_TGeoBBox_v), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("TGeoBBox_v",943,G__testCint_191_0_2, 105, G__get_linked_tagnum(&G__testCintLN_TGeoBBox_v), -1, 0, 4, 1, 1, 0, 
"d - 'Double_t' 0 - dx d - 'Double_t' 0 - dy "
"d - 'Double_t' 0 - dz D - 'Double_t' 0 '0' origin", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("TGeoBBox_v",943,G__testCint_191_0_3, 105, G__get_linked_tagnum(&G__testCintLN_TGeoBBox_v), -1, 0, 5, 1, 1, 0, 
"C - - 10 - name d - 'Double_t' 0 - dx "
"d - 'Double_t' 0 - dy d - 'Double_t' 0 - dz "
"D - 'Double_t' 0 '0' origin", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("TGeoBBox_v",943,G__testCint_191_0_4, 105, G__get_linked_tagnum(&G__testCintLN_TGeoBBox_v), -1, 0, 1, 1, 1, 0, "D - 'Double_t' 0 - param", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("AreOverlapping",1439,G__testCint_191_0_5, 103, -1, G__defined_typename("Bool_t"), 0, 4, 3, 1, 0, 
"U 'TGeoBBox_v' - 10 - box1 U 'TGeoMatrix' - 10 - mat1 "
"U 'TGeoBBox_v' - 10 - box2 U 'TGeoMatrix' - 10 - mat2", (char*)NULL, (void*) G__func2void( (Bool_t (*)(const TGeoBBox_v*, const TGeoMatrix*, const TGeoBBox_v*, const TGeoMatrix*))(&TGeoBBox_v::AreOverlapping) ), 0);
   G__memfunc_setup("ComputeNormal",1350,(G__InterfaceMethod) NULL,121, -1, -1, 0, 3, 1, 1, 0, 
"D - 'Double_t' 0 - point D - 'Double_t' 0 - dir "
"D - 'Double_t' 0 - norm", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Contains",831,(G__InterfaceMethod) NULL,103, -1, G__defined_typename("Bool_t"), 0, 1, 1, 1, 8, "D - 'Double_t' 0 - point", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Contains_v",1044,G__testCint_191_0_8, 121, -1, -1, 0, 3, 1, 1, 8, 
"D - 'Double_t' 10 - point G - 'Bool_t' 0 - isin "
"i - 'Int_t' 0 - np", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Contains",831,G__testCint_191_0_9, 103, -1, G__defined_typename("Bool_t"), 0, 5, 3, 1, 0, 
"D - 'Double_t' 10 - point d - 'Double_t' 0 - dx "
"d - 'Double_t' 0 - dy d - 'Double_t' 0 - dz "
"D - 'Double_t' 10 - origin", (char*)NULL, (void*) G__func2void( (Bool_t (*)(const Double_t*, Double_t, Double_t, Double_t, const Double_t*))(&TGeoBBox_v::Contains) ), 0);
   G__memfunc_setup("CouldBeCrossed",1393,(G__InterfaceMethod) NULL,103, -1, G__defined_typename("Bool_t"), 0, 2, 1, 1, 8, 
"D - 'Double_t' 0 - point D - 'Double_t' 0 - dir", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("DistancetoPrimitive",1991,(G__InterfaceMethod) NULL,105, -1, G__defined_typename("Int_t"), 0, 2, 1, 1, 0, 
"i - 'Int_t' 0 - px i - 'Int_t' 0 - py", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("DistFromInside",1412,(G__InterfaceMethod) NULL,100, -1, G__defined_typename("Double_t"), 0, 5, 1, 1, 8, 
"D - 'Double_t' 0 - point D - 'Double_t' 0 - dir "
"i - 'Int_t' 0 '1' iact d - 'Double_t' 0 'TGeoShape::Big()' step "
"D - 'Double_t' 0 '0' safe", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("DistFromInside",1412,G__testCint_191_0_13, 100, -1, G__defined_typename("Double_t"), 0, 7, 3, 1, 0, 
"D - 'Double_t' 10 - point D - 'Double_t' 10 - dir "
"d - 'Double_t' 0 - dx d - 'Double_t' 0 - dy "
"d - 'Double_t' 0 - dz D - 'Double_t' 10 - origin "
"d - 'Double_t' 0 'TGeoShape::Big()' stepmax", (char*)NULL, (void*) G__func2void( (Double_t (*)(const Double_t*, const Double_t*, Double_t, Double_t, Double_t, const Double_t*, Double_t))(&TGeoBBox_v::DistFromInside) ), 0);
   G__memfunc_setup("DistFromOutside",1541,(G__InterfaceMethod) NULL,100, -1, G__defined_typename("Double_t"), 0, 5, 1, 1, 8, 
"D - 'Double_t' 0 - point D - 'Double_t' 0 - dir "
"i - 'Int_t' 0 '1' iact d - 'Double_t' 0 'TGeoShape::Big()' step "
"D - 'Double_t' 0 '0' safe", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("DistFromOutside",1541,G__testCint_191_0_15, 100, -1, G__defined_typename("Double_t"), 0, 7, 3, 1, 0, 
"D - 'Double_t' 10 - point D - 'Double_t' 10 - dir "
"d - 'Double_t' 0 - dx d - 'Double_t' 0 - dy "
"d - 'Double_t' 0 - dz D - 'Double_t' 10 - origin "
"d - 'Double_t' 0 'TGeoShape::Big()' stepmax", (char*)NULL, (void*) G__func2void( (Double_t (*)(const Double_t*, const Double_t*, Double_t, Double_t, Double_t, const Double_t*, Double_t))(&TGeoBBox_v::DistFromOutside) ), 0);
   G__memfunc_setup("GetPointsOnFacet",1597,(G__InterfaceMethod) NULL,103, -1, G__defined_typename("Bool_t"), 0, 3, 1, 1, 8, 
"i - 'Int_t' 0 - index i - 'Int_t' 0 - npoints "
"D - 'Double_t' 0 - array", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("GetPointsOnSegments",1952,(G__InterfaceMethod) NULL,103, -1, G__defined_typename("Bool_t"), 0, 2, 1, 1, 8, 
"i - 'Int_t' 0 - npoints D - 'Double_t' 0 - array", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Safety",620,(G__InterfaceMethod) NULL,100, -1, G__defined_typename("Double_t"), 0, 2, 1, 1, 8, 
"D - 'Double_t' 0 - point g - 'Bool_t' 0 'kTRUE' in", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Class",502,G__testCint_191_0_19, 85, G__get_linked_tagnum(&G__testCintLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&TGeoBBox_v::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__testCint_191_0_20, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&TGeoBBox_v::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__testCint_191_0_21, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&TGeoBBox_v::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__testCint_191_0_22, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&TGeoBBox_v::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,(G__InterfaceMethod) NULL,85, G__get_linked_tagnum(&G__testCintLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TMemberInspector' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__testCint_191_0_26, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - ClassDef_StreamerNVirtual_b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__testCint_191_0_27, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&TGeoBBox_v::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__testCint_191_0_28, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&TGeoBBox_v::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__testCint_191_0_29, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&TGeoBBox_v::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__testCint_191_0_30, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&TGeoBBox_v::DeclFileLine) ), 0);
   // automatic copy constructor
   G__memfunc_setup("TGeoBBox_v", 943, G__testCint_191_0_31, (int) ('i'), G__get_linked_tagnum(&G__testCintLN_TGeoBBox_v), -1, 0, 1, 1, 1, 0, "u 'TGeoBBox_v' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~TGeoBBox_v", 1069, G__testCint_191_0_32, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 1);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__testCint_191_0_33, (int) ('u'), G__get_linked_tagnum(&G__testCintLN_TGeoBBox_v), -1, 1, 1, 1, 1, 0, "u 'TGeoBBox_v' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfunctestCint() {
}

/*********************************************************
* Global variable information setup for each class
*********************************************************/
static void G__cpp_setup_global0() {

   /* Setting up global variables */
   G__resetplocal();

}

static void G__cpp_setup_global1() {
}

static void G__cpp_setup_global2() {

   G__resetglobalenv();
}
extern "C" void G__cpp_setup_globaltestCint() {
  G__cpp_setup_global0();
  G__cpp_setup_global1();
  G__cpp_setup_global2();
}

/*********************************************************
* Global function information setup for each class
*********************************************************/
static void G__cpp_setup_func0() {
   G__lastifuncposition();

}

static void G__cpp_setup_func1() {
}

static void G__cpp_setup_func2() {
}

static void G__cpp_setup_func3() {
}

static void G__cpp_setup_func4() {
}

static void G__cpp_setup_func5() {
}

static void G__cpp_setup_func6() {
}

static void G__cpp_setup_func7() {
}

static void G__cpp_setup_func8() {
}

static void G__cpp_setup_func9() {
}

static void G__cpp_setup_func10() {
}

static void G__cpp_setup_func11() {
}

static void G__cpp_setup_func12() {
}

static void G__cpp_setup_func13() {

   G__resetifuncposition();
}

extern "C" void G__cpp_setup_functestCint() {
  G__cpp_setup_func0();
  G__cpp_setup_func1();
  G__cpp_setup_func2();
  G__cpp_setup_func3();
  G__cpp_setup_func4();
  G__cpp_setup_func5();
  G__cpp_setup_func6();
  G__cpp_setup_func7();
  G__cpp_setup_func8();
  G__cpp_setup_func9();
  G__cpp_setup_func10();
  G__cpp_setup_func11();
  G__cpp_setup_func12();
  G__cpp_setup_func13();
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__testCintLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__testCintLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__testCintLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__testCintLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__testCintLN_TNamed = { "TNamed" , 99 , -1 };
G__linked_taginfo G__testCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__testCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__testCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__testCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__testCintLN_TGeoBBox = { "TGeoBBox" , 99 , -1 };
G__linked_taginfo G__testCintLN_TGeoMatrix = { "TGeoMatrix" , 99 , -1 };
G__linked_taginfo G__testCintLN_TGeoShape = { "TGeoShape" , 99 , -1 };
G__linked_taginfo G__testCintLN_TGeoBBox_v = { "TGeoBBox_v" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtabletestCint() {
  G__testCintLN_TClass.tagnum = -1 ;
  G__testCintLN_TBuffer.tagnum = -1 ;
  G__testCintLN_TMemberInspector.tagnum = -1 ;
  G__testCintLN_TObject.tagnum = -1 ;
  G__testCintLN_TNamed.tagnum = -1 ;
  G__testCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__testCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__testCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__testCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__testCintLN_TGeoBBox.tagnum = -1 ;
  G__testCintLN_TGeoMatrix.tagnum = -1 ;
  G__testCintLN_TGeoShape.tagnum = -1 ;
  G__testCintLN_TGeoBBox_v.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtabletestCint() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__testCintLN_TClass);
   G__get_linked_tagnum_fwd(&G__testCintLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__testCintLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__testCintLN_TObject);
   G__get_linked_tagnum_fwd(&G__testCintLN_TNamed);
   G__get_linked_tagnum_fwd(&G__testCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__testCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__testCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__testCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__testCintLN_TGeoBBox);
   G__get_linked_tagnum_fwd(&G__testCintLN_TGeoMatrix);
   G__get_linked_tagnum_fwd(&G__testCintLN_TGeoShape);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__testCintLN_TGeoBBox_v),sizeof(TGeoBBox_v),-1,324864,"box primitive",G__setup_memvarTGeoBBox_v,G__setup_memfuncTGeoBBox_v);
}
extern "C" void G__cpp_setuptestCint(void) {
  G__check_setup_version(30051515,"G__cpp_setuptestCint()");
  G__set_cpp_environmenttestCint();
  G__cpp_setup_tagtabletestCint();

  G__cpp_setup_inheritancetestCint();

  G__cpp_setup_typetabletestCint();

  G__cpp_setup_memvartestCint();

  G__cpp_setup_memfunctestCint();
  G__cpp_setup_globaltestCint();
  G__cpp_setup_functestCint();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfunctestCint();
  return;
}
class G__cpp_setup_inittestCint {
  public:
    G__cpp_setup_inittestCint() { G__add_setup_func("testCint",(G__incsetup)(&G__cpp_setuptestCint)); G__call_setup_funcs(); }
   ~G__cpp_setup_inittestCint() { G__remove_setup_func("testCint"); }
};
G__cpp_setup_inittestCint G__cpp_setup_initializertestCint;
