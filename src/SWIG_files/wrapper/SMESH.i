/*
Copyright 2008-2015 Thomas Paviot (tpaviot@gmail.com)


This file is part of pythonOCC.
pythonOCC is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

pythonOCC is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with pythonOCC.  If not, see <http://www.gnu.org/licenses/>.

*/
%module (package="OCC") SMESH

#pragma SWIG nowarn=504,325,503

%{
#ifdef WNT
#pragma warning(disable : 4716)
#endif
%}

%include ../common/CommonIncludes.i
%include ../common/ExceptionCatcher.i
%include ../common/FunctionTransformers.i
%include ../common/Operators.i


%include SMESH_headers.i


%pythoncode {
def register_handle(handle, base_object):
    """
    Inserts the handle into the base object to
    prevent memory corruption in certain cases
    """
    try:
        if base_object.IsKind("Standard_Transient"):
            base_object.thisHandle = handle
            base_object.thisown = False
    except:
        pass
};

/* typedefs */
typedef std::map <const SMDS_MeshElement * , std::list <const SMDS_MeshElement *>> TElemOfElemListMap;
typedef const SMDS_MeshElement * SMDS_MeshElementPtr;
typedef std::map <SMESH_TLink , const SMDS_MeshNode *>::iterator ItTLinkNode;
typedef boost::shared_ptr <NumericalFunctor> SMESH::Controls::NumericalFunctorPtr;
typedef NCollection_IndexedMap <TopoDS_Shape> SMESH_IndexedMapOfShape;
typedef boost::shared_ptr <SMESH_ComputeError> SMESH_ComputeErrorPtr;
typedef boost::shared_ptr <GroupColor> SMESH::Controls::GroupColorPtr;
typedef boost::shared_ptr <SMESH_OctreeNodeIterator> SMESH_OctreeNodeIteratorPtr;
typedef NCollection_Sequence <SMDS_MeshElementPtr> SMESH_SequenceOfElemPtr;
typedef void ( * PVF ) ( );
typedef SMDS_Iterator <SMESH_OctreeNode *> SMESH_OctreeNodeIterator;
typedef boost::shared_ptr <Predicate> SMESH::Controls::PredicatePtr;
typedef boost::shared_ptr<SMDS_Iterator <SMESH_subMesh *>> SMESH_subMeshIteratorPtr;
typedef boost::shared_ptr <Length2D> SMESH::Controls::Length2DPtr;
typedef boost::shared_ptr <ManifoldPart> SMESH::Controls::ManifoldPartPtr;
typedef std::map <SMESH_TLink , const SMDS_MeshNode *> TLinkNodeMap;
typedef boost::shared_ptr <FreeEdges> SMESH::Controls::FreeEdgesPtr;
typedef NCollection_IndexedDataMap <SMESH_IndexedMapOfShape , TopoDS_Shape> SMESH_IndexedDataMapOfShapeIndexedMapOfShape;
typedef boost::shared_ptr <LogicalNOT> SMESH::Controls::LogicalNOTPtr;
typedef boost::shared_ptr <LogicalBinary> SMESH::Controls::LogicalBinaryPtr;
typedef std::map<double , TNodeColumn> TParam2ColumnMap;
typedef SMESH_Hypothesis::Hypothesis_Status TAlgoStateErrorName;
typedef boost::shared_ptr <RangeOfIds> SMESH::Controls::RangeOfIdsPtr;
typedef SMESH_subMeshEventListener EventListener;
typedef pair<const SMDS_MeshNode * , const SMDS_MeshNode *> NLink;
typedef const SMDS_MeshNode * SMDS_MeshNodePtr;
typedef NCollection_Sequence <SMDS_MeshNodePtr> SMESH_SequenceOfNode;
typedef std::map<SMESH_subMesh * , std::vector <int>>::iterator MapShapeNbElemsItr;
typedef boost::shared_ptr <MultiConnection2D> SMESH::Controls::MultiConnection2DPtr;
typedef boost::shared_ptr <Comparator> SMESH::Controls::ComparatorPtr;
typedef std::map<SMESH_subMesh * , std::vector <int>> MapShapeNbElems;
typedef boost::shared_ptr <Functor> SMESH::Controls::FunctorPtr;
typedef SMESH_subMeshEventListenerData EventListenerData;
typedef boost::shared_ptr <ElemGeomType> SMESH::Controls::ElemGeomTypePtr;
typedef boost::shared_ptr <ElementsOnSurface> SMESH::Controls::ElementsOnSurfacePtr;
typedef boost::shared_ptr <LinearOrQuadratic> SMESH::Controls::LinearOrQuadraticPtr;
typedef std::map <const SMDS_MeshNode * , const SMDS_MeshNode *> TNodeNodeMap;
typedef std::set <int> TSetOfInt;
typedef boost::shared_ptr <EqualTo> SMESH::Controls::EqualToPtr;
typedef std::set<const SMDS_MeshElement * , TIDCompare> TIDSortedElemSet;
typedef std::vector <const SMDS_MeshNode *> TNodeColumn;
typedef boost::shared_ptr <ElementsOnShape> SMESH::Controls::ElementsOnShapePtr;
/* end typedefs declaration */

/* public enums */
enum SMESH_ComputeErrorName {
	COMPERR_OK = - 1,
	COMPERR_BAD_INPUT_MESH = - 2,
	COMPERR_STD_EXCEPTION = - 3,
	COMPERR_OCC_EXCEPTION = - 4,
	COMPERR_SLM_EXCEPTION = - 5,
	COMPERR_EXCEPTION = - 6,
	COMPERR_MEMORY_PB = - 7,
	COMPERR_ALGO_FAILED = - 8,
	COMPERR_BAD_SHAPE = - 9,
};

enum MeshDimension {
	MeshDim_0D = 0,
	MeshDim_1D = 1,
	MeshDim_2D = 2,
	MeshDim_3D = 3,
};

/* end public enums declaration */

%nodefaultctor SMESH_Block;
class SMESH_Block : public math_FunctionSetWithDerivatives {
/* public enums */
enum TShapeID {
	ID_NONE = 0,
	ID_V000 = 1,
	ID_V100 = 2,
	ID_V010 = 3,
	ID_V110 = 4,
	ID_V001 = 5,
	ID_V101 = 6,
	ID_V011 = 7,
	ID_V111 = 8,
	ID_Ex00 = 9,
	ID_Ex10 = 10,
	ID_Ex01 = 11,
	ID_Ex11 = 12,
	ID_E0y0 = 13,
	ID_E1y0 = 14,
	ID_E0y1 = 15,
	ID_E1y1 = 16,
	ID_E00z = 17,
	ID_E10z = 18,
	ID_E01z = 19,
	ID_E11z = 20,
	ID_Fxy0 = 21,
	ID_Fxy1 = 22,
	ID_Fx0z = 23,
	ID_Fx1z = 24,
	ID_F0yz = 25,
	ID_F1yz = 26,
	ID_Shell = 27,
};

enum  {
	ID_FirstV = ID_V000,
	ID_FirstE = ID_Ex00,
	ID_FirstF = ID_Fxy0,
};

/* end public enums declaration */

	public:
		%feature("compactdefaultargs") NbVertices;
		%feature("autodoc", "	:rtype: int
") NbVertices;
		static int NbVertices ();
		%feature("compactdefaultargs") NbEdges;
		%feature("autodoc", "	:rtype: int
") NbEdges;
		static int NbEdges ();
		%feature("compactdefaultargs") NbFaces;
		%feature("autodoc", "	:rtype: int
") NbFaces;
		static int NbFaces ();
		%feature("compactdefaultargs") NbSubShapes;
		%feature("autodoc", "	:rtype: int
") NbSubShapes;
		static int NbSubShapes ();
		%feature("compactdefaultargs") IsVertexID;
		%feature("autodoc", "	:param theShapeID:
	:type theShapeID: int
	:rtype: inline bool
") IsVertexID;
		static inline bool IsVertexID (int theShapeID);
		%feature("compactdefaultargs") IsEdgeID;
		%feature("autodoc", "	:param theShapeID:
	:type theShapeID: int
	:rtype: inline bool
") IsEdgeID;
		static inline bool IsEdgeID (int theShapeID);
		%feature("compactdefaultargs") IsFaceID;
		%feature("autodoc", "	:param theShapeID:
	:type theShapeID: int
	:rtype: inline bool
") IsFaceID;
		static inline bool IsFaceID (int theShapeID);
		%feature("compactdefaultargs") ShapeIndex;
		%feature("autodoc", "	:param theShapeID:
	:type theShapeID: int
	:rtype: int
") ShapeIndex;
		static int ShapeIndex (int theShapeID);
		%feature("compactdefaultargs") GetFaceEdgesIDs;
		%feature("autodoc", "	:param faceID:
	:type faceID: int
	:param edgeVec:
	:type edgeVec: std::vector< int> &
	:rtype: void
") GetFaceEdgesIDs;
		static void GetFaceEdgesIDs (const int faceID,std::vector< int> & edgeVec);
		%feature("compactdefaultargs") GetEdgeVertexIDs;
		%feature("autodoc", "	:param edgeID:
	:type edgeID: int
	:param vertexVec:
	:type vertexVec: std::vector< int> &
	:rtype: void
") GetEdgeVertexIDs;
		static void GetEdgeVertexIDs (const int edgeID,std::vector< int> & vertexVec);
		%feature("compactdefaultargs") GetCoordIndOnEdge;
		%feature("autodoc", "	:param theEdgeID:
	:type theEdgeID: int
	:rtype: int
") GetCoordIndOnEdge;
		static int GetCoordIndOnEdge (const int theEdgeID);
		%feature("compactdefaultargs") GetShapeCoef;
		%feature("autodoc", "	:param theShapeID:
	:type theShapeID: int
	:rtype: double *
") GetShapeCoef;
		static double * GetShapeCoef (const int theShapeID);
		%feature("compactdefaultargs") GetShapeIDByParams;
		%feature("autodoc", "	:param theParams:
	:type theParams: gp_XYZ
	:rtype: int
") GetShapeIDByParams;
		static int GetShapeIDByParams (const gp_XYZ & theParams);
		%feature("compactdefaultargs") DumpShapeID;
		%feature("autodoc", "	:param theBlockShapeID:
	:type theBlockShapeID: int
	:param stream:
	:type stream: std::ostream &
	:rtype: std::ostream
") DumpShapeID;
		static std::ostream & DumpShapeID (const int theBlockShapeID,std::ostream & stream);
		%feature("compactdefaultargs") SMESH_Block;
		%feature("autodoc", "	:rtype: None
") SMESH_Block;
		 SMESH_Block ();
		%feature("compactdefaultargs") LoadBlockShapes;
		%feature("autodoc", "	:param theShell:
	:type theShell: TopoDS_Shell &
	:param theVertex000:
	:type theVertex000: TopoDS_Vertex &
	:param theVertex001:
	:type theVertex001: TopoDS_Vertex &
	:param theShapeIDMap:
	:type theShapeIDMap: TopTools_IndexedMapOfOrientedShape &
	:rtype: bool
") LoadBlockShapes;
		bool LoadBlockShapes (const TopoDS_Shell & theShell,const TopoDS_Vertex & theVertex000,const TopoDS_Vertex & theVertex001,TopTools_IndexedMapOfOrientedShape & theShapeIDMap);
		%feature("compactdefaultargs") LoadBlockShapes;
		%feature("autodoc", "	:param theShapeIDMap:
	:type theShapeIDMap: TopTools_IndexedMapOfOrientedShape &
	:rtype: bool
") LoadBlockShapes;
		bool LoadBlockShapes (const TopTools_IndexedMapOfOrientedShape & theShapeIDMap);
		%feature("compactdefaultargs") LoadMeshBlock;
		%feature("autodoc", "	:param theVolume:
	:type theVolume: SMDS_MeshVolume *
	:param theNode000Index:
	:type theNode000Index: int
	:param theNode001Index:
	:type theNode001Index: int
	:param theOrderedNodes:
	:type theOrderedNodes: std::vector< SMDS_MeshNode *> &
	:rtype: bool
") LoadMeshBlock;
		bool LoadMeshBlock (const SMDS_MeshVolume * theVolume,const int theNode000Index,const int theNode001Index,std::vector<const SMDS_MeshNode *> & theOrderedNodes);
		%feature("compactdefaultargs") LoadFace;
		%feature("autodoc", "	:param theFace:
	:type theFace: TopoDS_Face &
	:param theFaceID:
	:type theFaceID: int
	:param theShapeIDMap:
	:type theShapeIDMap: TopTools_IndexedMapOfOrientedShape &
	:rtype: bool
") LoadFace;
		bool LoadFace (const TopoDS_Face & theFace,const int theFaceID,const TopTools_IndexedMapOfOrientedShape & theShapeIDMap);
		%feature("compactdefaultargs") Insert;
		%feature("autodoc", "	:param theShape:
	:type theShape: TopoDS_Shape &
	:param theShapeID:
	:type theShapeID: int
	:param theShapeIDMap:
	:type theShapeIDMap: TopTools_IndexedMapOfOrientedShape &
	:rtype: bool
") Insert;
		static bool Insert (const TopoDS_Shape & theShape,const int theShapeID,TopTools_IndexedMapOfOrientedShape & theShapeIDMap);
		%feature("compactdefaultargs") FindBlockShapes;
		%feature("autodoc", "	:param theShell:
	:type theShell: TopoDS_Shell &
	:param theVertex000:
	:type theVertex000: TopoDS_Vertex &
	:param theVertex001:
	:type theVertex001: TopoDS_Vertex &
	:param theShapeIDMap:
	:type theShapeIDMap: TopTools_IndexedMapOfOrientedShape &
	:rtype: bool
") FindBlockShapes;
		static bool FindBlockShapes (const TopoDS_Shell & theShell,const TopoDS_Vertex & theVertex000,const TopoDS_Vertex & theVertex001,TopTools_IndexedMapOfOrientedShape & theShapeIDMap);
		%feature("compactdefaultargs") VertexPoint;
		%feature("autodoc", "	:param theVertexID:
	:type theVertexID: int
	:param thePoint:
	:type thePoint: gp_XYZ
	:rtype: bool
") VertexPoint;
		bool VertexPoint (const int theVertexID,gp_XYZ & thePoint);
		%feature("compactdefaultargs") EdgePoint;
		%feature("autodoc", "	:param theEdgeID:
	:type theEdgeID: int
	:param theParams:
	:type theParams: gp_XYZ
	:param thePoint:
	:type thePoint: gp_XYZ
	:rtype: bool
") EdgePoint;
		bool EdgePoint (const int theEdgeID,const gp_XYZ & theParams,gp_XYZ & thePoint);
		%feature("compactdefaultargs") EdgeU;
		%feature("autodoc", "	:param theEdgeID:
	:type theEdgeID: int
	:param theParams:
	:type theParams: gp_XYZ
	:param theU:
	:type theU: double &
	:rtype: bool
") EdgeU;
		bool EdgeU (const int theEdgeID,const gp_XYZ & theParams,Standard_Real &OutValue);
		%feature("compactdefaultargs") FacePoint;
		%feature("autodoc", "	:param theFaceID:
	:type theFaceID: int
	:param theParams:
	:type theParams: gp_XYZ
	:param thePoint:
	:type thePoint: gp_XYZ
	:rtype: bool
") FacePoint;
		bool FacePoint (const int theFaceID,const gp_XYZ & theParams,gp_XYZ & thePoint);
		%feature("compactdefaultargs") FaceUV;
		%feature("autodoc", "	:param theFaceID:
	:type theFaceID: int
	:param theParams:
	:type theParams: gp_XYZ
	:param theUV:
	:type theUV: gp_XY
	:rtype: bool
") FaceUV;
		bool FaceUV (const int theFaceID,const gp_XYZ & theParams,gp_XY & theUV);
		%feature("compactdefaultargs") ShellPoint;
		%feature("autodoc", "	:param theParams:
	:type theParams: gp_XYZ
	:param thePoint:
	:type thePoint: gp_XYZ
	:rtype: bool
") ShellPoint;
		bool ShellPoint (const gp_XYZ & theParams,gp_XYZ & thePoint);
		%feature("compactdefaultargs") ShellPoint;
		%feature("autodoc", "	:param theParams:
	:type theParams: gp_XYZ
	:param thePointOnShape:
	:type thePointOnShape: std::vector<gp_XYZ>
	:param thePoint:
	:type thePoint: gp_XYZ
	:rtype: bool
") ShellPoint;
		static bool ShellPoint (const gp_XYZ & theParams,const std::vector<gp_XYZ> & thePointOnShape,gp_XYZ & thePoint);
		%feature("compactdefaultargs") ComputeParameters;
		%feature("autodoc", "	:param thePoint:
	:type thePoint: gp_Pnt
	:param theParams:
	:type theParams: gp_XYZ
	:param theShapeID: default value is ID_Shell
	:type theShapeID: int
	:param theParamsHint: default value is gp_XYZ(-1,-1,-1)
	:type theParamsHint: gp_XYZ
	:rtype: bool
") ComputeParameters;
		bool ComputeParameters (const gp_Pnt & thePoint,gp_XYZ & theParams,const int theShapeID = ID_Shell,const gp_XYZ & theParamsHint = gp_XYZ(-1,-1,-1));
		%feature("compactdefaultargs") VertexParameters;
		%feature("autodoc", "	:param theVertexID:
	:type theVertexID: int
	:param theParams:
	:type theParams: gp_XYZ
	:rtype: bool
") VertexParameters;
		bool VertexParameters (const int theVertexID,gp_XYZ & theParams);
		%feature("compactdefaultargs") EdgeParameters;
		%feature("autodoc", "	:param theEdgeID:
	:type theEdgeID: int
	:param theU:
	:type theU: double
	:param theParams:
	:type theParams: gp_XYZ
	:rtype: bool
") EdgeParameters;
		bool EdgeParameters (const int theEdgeID,const double theU,gp_XYZ & theParams);
		%feature("compactdefaultargs") IsForwardEdge;
		%feature("autodoc", "	:param theEdge:
	:type theEdge: TopoDS_Edge &
	:param theShapeIDMap:
	:type theShapeIDMap: TopTools_IndexedMapOfOrientedShape &
	:rtype: bool
") IsForwardEdge;
		static bool IsForwardEdge (const TopoDS_Edge & theEdge,const TopTools_IndexedMapOfOrientedShape & theShapeIDMap);
		%feature("compactdefaultargs") GetOrderedEdges;
		%feature("autodoc", "	:param theFace:
	:type theFace: TopoDS_Face &
	:param theFirstVertex:
	:type theFirstVertex: TopoDS_Vertex
	:param theEdges:
	:type theEdges: std::list< TopoDS_Edge> &
	:param theNbVertexInWires:
	:type theNbVertexInWires: std::list< int> &
	:param theShapeAnalysisAlgo: default value is false
	:type theShapeAnalysisAlgo: bool
	:rtype: int
") GetOrderedEdges;
		static int GetOrderedEdges (const TopoDS_Face & theFace,TopoDS_Vertex theFirstVertex,std::list< TopoDS_Edge> & theEdges,std::list< int> & theNbVertexInWires,const bool theShapeAnalysisAlgo = false);
		%feature("compactdefaultargs") NbVariables;
		%feature("autodoc", "	:rtype: int
") NbVariables;
		Standard_Integer NbVariables ();
		%feature("compactdefaultargs") NbEquations;
		%feature("autodoc", "	:rtype: int
") NbEquations;
		Standard_Integer NbEquations ();
		%feature("compactdefaultargs") Value;
		%feature("autodoc", "	:param X:
	:type X: math_Vector &
	:param F:
	:type F: math_Vector &
	:rtype: bool
") Value;
		Standard_Boolean Value (const math_Vector & X,math_Vector & F);
		%feature("compactdefaultargs") Derivatives;
		%feature("autodoc", "	:param X:
	:type X: math_Vector &
	:param D:
	:type D: math_Matrix &
	:rtype: bool
") Derivatives;
		Standard_Boolean Derivatives (const math_Vector & X,math_Matrix & D);
		%feature("compactdefaultargs") Values;
		%feature("autodoc", "	:param X:
	:type X: math_Vector &
	:param F:
	:type F: math_Vector &
	:param D:
	:type D: math_Matrix &
	:rtype: bool
") Values;
		Standard_Boolean Values (const math_Vector & X,math_Vector & F,math_Matrix & D);
		%feature("compactdefaultargs") GetStateNumber;
		%feature("autodoc", "	:rtype: int
") GetStateNumber;
		Standard_Integer GetStateNumber ();
};


%nodefaultctor SMESH_ComputeError;
class SMESH_ComputeError {
	public:
		%feature("compactdefaultargs") New;
		%feature("autodoc", "	* //!< to explain COMPERR_BAD_INPUT_MESH

	:param error: default value is COMPERR_OK
	:type error: int
	:param comment: default value is ""
	:type comment: std::string
	:param algo: default value is 0
	:type algo: SMESH_Algo *
	:rtype: SMESH_ComputeErrorPtr
") New;
		static SMESH_ComputeErrorPtr New (int error = COMPERR_OK,std::string comment = "",const SMESH_Algo * algo = 0);
		%feature("compactdefaultargs") SMESH_ComputeError;
		%feature("autodoc", "	:param error: default value is COMPERR_OK
	:type error: int
	:param comment: default value is ""
	:type comment: std::string
	:param algo: default value is 0
	:type algo: SMESH_Algo *
	:rtype: None
") SMESH_ComputeError;
		 SMESH_ComputeError (int error = COMPERR_OK,std::string comment = "",const SMESH_Algo * algo = 0);
		%feature("compactdefaultargs") IsOK;
		%feature("autodoc", "	:rtype: bool
") IsOK;
		bool IsOK ();
		%feature("compactdefaultargs") IsCommon;
		%feature("autodoc", "	:rtype: bool
") IsCommon;
		bool IsCommon ();
		%feature("compactdefaultargs") CommonName;
		%feature("autodoc", "	:rtype: inline std::string
") CommonName;
		inline std::string CommonName ();
};


%nodefaultctor SMESH_ElementSearcher;
class SMESH_ElementSearcher {
	public:
		%feature("compactdefaultargs") FindElementsByPoint;
		%feature("autodoc", "	:param point:
	:type point: gp_Pnt
	:param type:
	:type type: SMDSAbs_ElementType
	:param foundElems:
	:type foundElems: std::vector<  SMDS_MeshElement *> &
	:rtype: None
") FindElementsByPoint;
		void FindElementsByPoint (const gp_Pnt & point,SMDSAbs_ElementType type,std::vector< const SMDS_MeshElement *> & foundElems);
};


%nodefaultctor SMESH_Exception;
class SMESH_Exception : public std::exception {
	public:
		%feature("compactdefaultargs") SMESH_Exception;
		%feature("autodoc", "	:param text:
	:type text: char *
	:param fileName: default value is 0
	:type fileName: char *
	:param lineNumber: default value is 0
	:type lineNumber: unsigned int
	:rtype: None
") SMESH_Exception;
		 SMESH_Exception (const char * text,const char * fileName = 0,const unsigned int lineNumber = 0);
		%feature("compactdefaultargs") SMESH_Exception;
		%feature("autodoc", "	:param ex:
	:type ex: SMESH_Exception &
	:rtype: None
") SMESH_Exception;
		 SMESH_Exception (const SMESH_Exception & ex);
		%feature("compactdefaultargs") what;
		%feature("autodoc", "	:param :
	:type : void
	:rtype: char *
") what;
		const char * what (void );
};


%nodefaultctor SMESH_Gen;
class SMESH_Gen {
	public:
		%feature("compactdefaultargs") SMESH_Gen;
		%feature("autodoc", "	:rtype: None
") SMESH_Gen;
		 SMESH_Gen ();
		%feature("compactdefaultargs") CreateMesh;
		%feature("autodoc", "	:param theStudyId:
	:type theStudyId: int
	:param theIsEmbeddedMode:
	:type theIsEmbeddedMode: bool
	:rtype: SMESH_Mesh *
") CreateMesh;
		SMESH_Mesh * CreateMesh (int theStudyId,bool theIsEmbeddedMode);
		%feature("compactdefaultargs") Compute;
		%feature("autodoc", "	* /*! * \brief Computes aMesh on aShape * \param anUpward - compute from vertices up to more complex shape (internal usage) * \param aDim - upper level dimension of the mesh computation * \param aShapesId - list of shapes with computed mesh entities (elements or nodes) * etval bool - true if none submesh failed to compute */

	:param aMesh:
	:type aMesh: SMESH_Mesh &
	:param aShape:
	:type aShape: TopoDS_Shape &
	:param anUpward: default value is false
	:type anUpward: bool
	:param aDim: default value is MeshDim_3D
	:type aDim: MeshDimension
	:param aShapesId: default value is 0
	:type aShapesId: TSetOfInt *
	:rtype: bool
") Compute;
		bool Compute (SMESH_Mesh & aMesh,const TopoDS_Shape & aShape,const bool anUpward = false,const MeshDimension aDim = MeshDim_3D,TSetOfInt * aShapesId = 0);
		%feature("compactdefaultargs") Evaluate;
		%feature("autodoc", "	* /*! * \brief evaluates size of prospective mesh on a shape * \param aMesh - the mesh * \param aShape - the shape * \param aResMap - map for prospective numbers of elements * etval bool - is a success */

	:param aMesh:
	:type aMesh: SMESH_Mesh &
	:param aShape:
	:type aShape: TopoDS_Shape &
	:param aResMap:
	:type aResMap: MapShapeNbElems &
	:param anUpward: default value is false
	:type anUpward: bool
	:param aShapesId: default value is 0
	:type aShapesId: TSetOfInt *
	:rtype: bool
") Evaluate;
		bool Evaluate (SMESH_Mesh & aMesh,const TopoDS_Shape & aShape,MapShapeNbElems & aResMap,const bool anUpward = false,TSetOfInt * aShapesId = 0);
		%feature("compactdefaultargs") CheckAlgoState;
		%feature("autodoc", "	:param aMesh:
	:type aMesh: SMESH_Mesh &
	:param aShape:
	:type aShape: TopoDS_Shape &
	:rtype: bool
") CheckAlgoState;
		bool CheckAlgoState (SMESH_Mesh & aMesh,const TopoDS_Shape & aShape);
		%feature("compactdefaultargs") SetBoundaryBoxSegmentation;
		%feature("autodoc", "	* /*! * \brief Sets number of segments per diagonal of boundary box of geometry by which * default segment length of appropriate 1D hypotheses is defined */

	:param theNbSegments:
	:type theNbSegments: int
	:rtype: None
") SetBoundaryBoxSegmentation;
		void SetBoundaryBoxSegmentation (int theNbSegments);
		%feature("compactdefaultargs") GetBoundaryBoxSegmentation;
		%feature("autodoc", "	:rtype: int
") GetBoundaryBoxSegmentation;
		int GetBoundaryBoxSegmentation ();
		%feature("compactdefaultargs") SetDefaultNbSegments;
		%feature("autodoc", "	* /*! * \brief Sets default number of segments per edge */

	:param nb:
	:type nb: int
	:rtype: None
") SetDefaultNbSegments;
		void SetDefaultNbSegments (int nb);
		%feature("compactdefaultargs") GetDefaultNbSegments;
		%feature("autodoc", "	:rtype: int
") GetDefaultNbSegments;
		int GetDefaultNbSegments ();
		%feature("compactdefaultargs") GetAlgoState;
		%feature("autodoc", "	:param aMesh:
	:type aMesh: SMESH_Mesh &
	:param aShape:
	:type aShape: TopoDS_Shape &
	:param theErrors:
	:type theErrors: std::list< SMESH_Gen::TAlgoStateError> &
	:rtype: bool
") GetAlgoState;
		bool GetAlgoState (SMESH_Mesh & aMesh,const TopoDS_Shape & aShape,std::list< SMESH_Gen::TAlgoStateError> & theErrors);
		%feature("compactdefaultargs") GetStudyContext;
		%feature("autodoc", "	:param studyId:
	:type studyId: int
	:rtype: StudyContextStruct *
") GetStudyContext;
		StudyContextStruct * GetStudyContext (int studyId);
		%feature("compactdefaultargs") GetShapeDim;
		%feature("autodoc", "	:param aShapeType:
	:type aShapeType: TopAbs_ShapeEnum &
	:rtype: int
") GetShapeDim;
		static int GetShapeDim (const TopAbs_ShapeEnum & aShapeType);
		%feature("compactdefaultargs") GetShapeDim;
		%feature("autodoc", "	:param aShape:
	:type aShape: TopoDS_Shape &
	:rtype: int
") GetShapeDim;
		static int GetShapeDim (const TopoDS_Shape & aShape);
		%feature("compactdefaultargs") GetAlgo;
		%feature("autodoc", "	:param aMesh:
	:type aMesh: SMESH_Mesh &
	:param aShape:
	:type aShape: TopoDS_Shape &
	:param assignedTo: default value is 0
	:type assignedTo: TopoDS_Shape *
	:rtype: SMESH_Algo *
") GetAlgo;
		SMESH_Algo * GetAlgo (SMESH_Mesh & aMesh,const TopoDS_Shape & aShape,TopoDS_Shape * assignedTo = 0);
		%feature("compactdefaultargs") IsGlobalHypothesis;
		%feature("autodoc", "	:param theHyp:
	:type theHyp: SMESH_Hypothesis *
	:param aMesh:
	:type aMesh: SMESH_Mesh &
	:rtype: bool
") IsGlobalHypothesis;
		static bool IsGlobalHypothesis (const SMESH_Hypothesis * theHyp,SMESH_Mesh & aMesh);
		%feature("compactdefaultargs") GetANewId;
		%feature("autodoc", "	:rtype: int
") GetANewId;
		int GetANewId ();
};


%nodefaultctor SMESH_Group;
class SMESH_Group {
	public:
		%feature("compactdefaultargs") SMESH_Group;
		%feature("autodoc", "	:param theID:
	:type theID: int
	:param theMesh:
	:type theMesh: SMESH_Mesh *
	:param theType:
	:type theType: SMDSAbs_ElementType
	:param theName:
	:type theName: char *
	:param theShape: default value is TopoDS_Shape()
	:type theShape: TopoDS_Shape &
	:rtype: None
") SMESH_Group;
		 SMESH_Group (int theID,const SMESH_Mesh * theMesh,const SMDSAbs_ElementType theType,const char * theName,const TopoDS_Shape & theShape = TopoDS_Shape());
		%feature("compactdefaultargs") SetName;
		%feature("autodoc", "	:param theName:
	:type theName: char *
	:rtype: None
") SetName;
		void SetName (const char * theName);
		%feature("compactdefaultargs") GetName;
		%feature("autodoc", "	:rtype: char *
") GetName;
		const char * GetName ();
		%feature("compactdefaultargs") GetGroupDS;
		%feature("autodoc", "	:rtype: SMESHDS_GroupBase *
") GetGroupDS;
		SMESHDS_GroupBase * GetGroupDS ();
};


%nodefaultctor SMESH_HypoPredicate;
class SMESH_HypoPredicate {
	public:
		%feature("compactdefaultargs") IsOk;
		%feature("autodoc", "	:param aHyp:
	:type aHyp: SMESH_Hypothesis *
	:param aShape:
	:type aShape: TopoDS_Shape &
	:rtype: bool
") IsOk;
		bool IsOk (const SMESH_Hypothesis * aHyp,const TopoDS_Shape & aShape);
};


%nodefaultctor SMESH_Hypothesis;
class SMESH_Hypothesis : public SMESHDS_Hypothesis {
/* public enums */
enum Hypothesis_Status {
	HYP_OK = 0,
	HYP_MISSING = 1,
	HYP_CONCURENT = 2,
	HYP_BAD_PARAMETER = 3,
	HYP_HIDDEN_ALGO = 4,
	HYP_HIDING_ALGO = 5,
	HYP_UNKNOWN_FATAL = 6,
	HYP_INCOMPATIBLE = 7,
	HYP_NOTCONFORM = 8,
	HYP_ALREADY_EXIST = 9,
	HYP_BAD_DIM = 10,
	HYP_BAD_SUBSHAPE = 11,
	HYP_BAD_GEOMETRY = 12,
	HYP_NEED_SHAPE = 13,
};

/* end public enums declaration */

	public:
		%feature("compactdefaultargs") GetDim;
		%feature("autodoc", "	:rtype: int
") GetDim;
		int GetDim ();
		%feature("compactdefaultargs") GetStudyId;
		%feature("autodoc", "	:rtype: int
") GetStudyId;
		int GetStudyId ();
		%feature("compactdefaultargs") NotifySubMeshesHypothesisModification;
		%feature("autodoc", "	:rtype: None
") NotifySubMeshesHypothesisModification;
		void NotifySubMeshesHypothesisModification ();
		%feature("compactdefaultargs") GetShapeType;
		%feature("autodoc", "	:rtype: int
") GetShapeType;
		int GetShapeType ();
		%feature("compactdefaultargs") GetLibName;
		%feature("autodoc", "	:rtype: char *
") GetLibName;
		const char * GetLibName ();
		%feature("compactdefaultargs") SetLibName;
		%feature("autodoc", "	:param theLibName:
	:type theLibName: char *
	:rtype: None
") SetLibName;
		void SetLibName (const char * theLibName);
		%feature("compactdefaultargs") SetParameters;
		%feature("autodoc", "	:param theParameters:
	:type theParameters: char *
	:rtype: None
") SetParameters;
		void SetParameters (const char * theParameters);
		%feature("compactdefaultargs") GetParameters;
		%feature("autodoc", "	:rtype: char *
") GetParameters;
		char * GetParameters ();
		%feature("compactdefaultargs") SetLastParameters;
		%feature("autodoc", "	:param theParameters:
	:type theParameters: char *
	:rtype: None
") SetLastParameters;
		void SetLastParameters (const char * theParameters);
		%feature("compactdefaultargs") GetLastParameters;
		%feature("autodoc", "	:rtype: char *
") GetLastParameters;
		char * GetLastParameters ();
		%feature("compactdefaultargs") ClearParameters;
		%feature("autodoc", "	:rtype: None
") ClearParameters;
		void ClearParameters ();
		%feature("compactdefaultargs") SetParametersByMesh;
		%feature("autodoc", "	* /*! * \brief Initialize my parameter values by the mesh built on the geometry * \param theMesh - the built mesh * \param theShape - the geometry of interest * etval bool - true if parameter values have been successfully defined */

	:param theMesh:
	:type theMesh: SMESH_Mesh *
	:param theShape:
	:type theShape: TopoDS_Shape &
	:rtype: bool
") SetParametersByMesh;
		bool SetParametersByMesh (const SMESH_Mesh * theMesh,const TopoDS_Shape & theShape);
		%feature("compactdefaultargs") SetParametersByDefaults;
		%feature("autodoc", "	* /*! * \brief Initialize my parameter values by default parameters. * etval bool - true if parameter values have been successfully defined */

	:param dflts:
	:type dflts: SMESH_0D_Algo::TDefaults &
	:param theMesh: default value is 0
	:type theMesh: SMESH_Mesh *
	:rtype: bool
") SetParametersByDefaults;
		bool SetParametersByDefaults (const SMESH_0D_Algo::TDefaults & dflts,const SMESH_Mesh * theMesh = 0);
		%feature("compactdefaultargs") IsAuxiliary;
		%feature("autodoc", "	* /*! * \brief Return true if me is an auxiliary hypothesis * etval bool - auxiliary or not * * An auxiliary hypothesis is optional, i.e. an algorithm * can work without it and another hypothesis of the same * dimention can be assigned to the shape */

	:rtype: bool
") IsAuxiliary;
		bool IsAuxiliary ();
};


%nodefaultctor SMESH_Mesh;
class SMESH_Mesh {
	public:
		%feature("compactdefaultargs") SMESH_Mesh;
		%feature("autodoc", "	:param theLocalId:
	:type theLocalId: int
	:param theStudyId:
	:type theStudyId: int
	:param theGen:
	:type theGen: SMESH_Gen *
	:param theIsEmbeddedMode:
	:type theIsEmbeddedMode: bool
	:param theDocument:
	:type theDocument: SMESHDS_Document *
	:rtype: None
") SMESH_Mesh;
		 SMESH_Mesh (int theLocalId,int theStudyId,SMESH_Gen * theGen,bool theIsEmbeddedMode,SMESHDS_Document * theDocument);
		%feature("compactdefaultargs") ShapeToMesh;
		%feature("autodoc", "	* /*! * \brief Set geometry to be meshed */

	:param aShape:
	:type aShape: TopoDS_Shape &
	:rtype: None
") ShapeToMesh;
		void ShapeToMesh (const TopoDS_Shape & aShape);
		%feature("compactdefaultargs") GetShapeToMesh;
		%feature("autodoc", "	* /*! * \brief Return geometry to be meshed. (It may be a PseudoShape()!) */

	:rtype: TopoDS_Shape
") GetShapeToMesh;
		TopoDS_Shape GetShapeToMesh ();
		%feature("compactdefaultargs") HasShapeToMesh;
		%feature("autodoc", "	* /*! * \brief Return true if there is a geometry to be meshed, not PseudoShape() */

	:rtype: bool
") HasShapeToMesh;
		bool HasShapeToMesh ();
		%feature("compactdefaultargs") GetShapeDiagonalSize;
		%feature("autodoc", "	* /*! * \brief Return diagonal size of bounding box of shape to mesh. */

	:rtype: double
") GetShapeDiagonalSize;
		double GetShapeDiagonalSize ();
		%feature("compactdefaultargs") GetShapeDiagonalSize;
		%feature("autodoc", "	* /*! * \brief Return diagonal size of bounding box of a shape. */

	:param aShape:
	:type aShape: TopoDS_Shape &
	:rtype: double
") GetShapeDiagonalSize;
		static double GetShapeDiagonalSize (const TopoDS_Shape & aShape);
		%feature("compactdefaultargs") PseudoShape;
		%feature("autodoc", "	* /*! * \brief Return a solid which is returned by GetShapeToMesh() if * a real geometry to be meshed was not set */

	:rtype: TopoDS_Solid
") PseudoShape;
		static const TopoDS_Solid  PseudoShape ();
		%feature("compactdefaultargs") Clear;
		%feature("autodoc", "	* /*! * \brief Remove all nodes and elements */

	:rtype: None
") Clear;
		void Clear ();
		%feature("compactdefaultargs") ClearSubMesh;
		%feature("autodoc", "	* /*! * \brief Remove all nodes and elements of indicated shape */

	:param theShapeId:
	:type theShapeId: int
	:rtype: None
") ClearSubMesh;
		void ClearSubMesh (const int theShapeId);
		%feature("compactdefaultargs") UNVToMesh;
		%feature("autodoc", "	:param theFileName:
	:type theFileName: char *
	:rtype: int
") UNVToMesh;
		int UNVToMesh (const char * theFileName);
		%feature("compactdefaultargs") MEDToMesh;
		%feature("autodoc", "	* /*! * consult DriverMED_R_SMESHDS_Mesh::ReadStatus for returned value */

	:param theFileName:
	:type theFileName: char *
	:param theMeshName:
	:type theMeshName: char *
	:rtype: int
") MEDToMesh;
		int MEDToMesh (const char * theFileName,const char * theMeshName);
		%feature("compactdefaultargs") STLToMesh;
		%feature("autodoc", "	:param theFileName:
	:type theFileName: char *
	:rtype: int
") STLToMesh;
		int STLToMesh (const char * theFileName);
		%feature("compactdefaultargs") DATToMesh;
		%feature("autodoc", "	:param theFileName:
	:type theFileName: char *
	:rtype: int
") DATToMesh;
		int DATToMesh (const char * theFileName);
		%feature("compactdefaultargs") AddHypothesis;
		%feature("autodoc", "	:param aSubShape:
	:type aSubShape: TopoDS_Shape &
	:param anHypId:
	:type anHypId: int
	:rtype: SMESH_Hypothesis::Hypothesis_Status
") AddHypothesis;
		SMESH_Hypothesis::Hypothesis_Status AddHypothesis (const TopoDS_Shape & aSubShape,int anHypId);
		%feature("compactdefaultargs") RemoveHypothesis;
		%feature("autodoc", "	:param aSubShape:
	:type aSubShape: TopoDS_Shape &
	:param anHypId:
	:type anHypId: int
	:rtype: SMESH_Hypothesis::Hypothesis_Status
") RemoveHypothesis;
		SMESH_Hypothesis::Hypothesis_Status RemoveHypothesis (const TopoDS_Shape & aSubShape,int anHypId);
		%feature("compactdefaultargs") GetHypothesisList;
		%feature("autodoc", "	:param aSubShape:
	:type aSubShape: TopoDS_Shape &
	:rtype: std::list< SMESHDS_Hypothesis *>
") GetHypothesisList;
		const std::list<const SMESHDS_Hypothesis *> & GetHypothesisList (const TopoDS_Shape & aSubShape);
		%feature("compactdefaultargs") GetHypothesis;
		%feature("autodoc", "	:param aSubShape:
	:type aSubShape: TopoDS_Shape &
	:param aFilter:
	:type aFilter: SMESH_HypoFilter &
	:param andAncestors:
	:type andAncestors: bool
	:param assignedTo: default value is 0
	:type assignedTo: TopoDS_Shape *
	:rtype: SMESH_Hypothesis *
") GetHypothesis;
		const SMESH_Hypothesis * GetHypothesis (const TopoDS_Shape & aSubShape,const SMESH_HypoFilter & aFilter,const bool andAncestors,TopoDS_Shape * assignedTo = 0);
		%feature("compactdefaultargs") GetHypotheses;
		%feature("autodoc", "	:param aSubShape:
	:type aSubShape: TopoDS_Shape &
	:param aFilter:
	:type aFilter: SMESH_HypoFilter &
	:param aHypList:
	:type aHypList: std::list< SMESHDS_Hypothesis *> &
	:param andAncestors:
	:type andAncestors: bool
	:rtype: int
") GetHypotheses;
		int GetHypotheses (const TopoDS_Shape & aSubShape,const SMESH_HypoFilter & aFilter,std::list<const SMESHDS_Hypothesis *> & aHypList,const bool andAncestors);
		%feature("compactdefaultargs") GetLog;
		%feature("autodoc", "	:rtype: std::list<SMESHDS_Command *>
") GetLog;
		const std::list<SMESHDS_Command *> & GetLog ();
		%feature("compactdefaultargs") ClearLog;
		%feature("autodoc", "	:rtype: None
") ClearLog;
		void ClearLog ();
		%feature("compactdefaultargs") GetId;
		%feature("autodoc", "	:rtype: int
") GetId;
		int GetId ();
		%feature("compactdefaultargs") GetMeshDS;
		%feature("autodoc", "	:rtype: SMESHDS_Mesh *
") GetMeshDS;
		SMESHDS_Mesh * GetMeshDS ();
		%feature("compactdefaultargs") GetGen;
		%feature("autodoc", "	:rtype: SMESH_Gen *
") GetGen;
		SMESH_Gen * GetGen ();
		%feature("compactdefaultargs") GetSubMesh;
		%feature("autodoc", "	:param aSubShape:
	:type aSubShape: TopoDS_Shape &
	:rtype: SMESH_subMesh *
") GetSubMesh;
		SMESH_subMesh * GetSubMesh (const TopoDS_Shape & aSubShape);
		%feature("compactdefaultargs") GetSubMeshContaining;
		%feature("autodoc", "	:param aSubShape:
	:type aSubShape: TopoDS_Shape &
	:rtype: SMESH_subMesh *
") GetSubMeshContaining;
		SMESH_subMesh * GetSubMeshContaining (const TopoDS_Shape & aSubShape);
		%feature("compactdefaultargs") GetSubMeshContaining;
		%feature("autodoc", "	:param aShapeID:
	:type aShapeID: int
	:rtype: SMESH_subMesh *
") GetSubMeshContaining;
		SMESH_subMesh * GetSubMeshContaining (const int aShapeID);
		%feature("compactdefaultargs") GetGroupSubMeshesContaining;
		%feature("autodoc", "	* /*! * \brief Return submeshes of groups containing the given subshape */

	:param shape:
	:type shape: TopoDS_Shape &
	:rtype: std::list<SMESH_subMesh *>
") GetGroupSubMeshesContaining;
		std::list<SMESH_subMesh *> GetGroupSubMeshesContaining (const TopoDS_Shape & shape);
		%feature("compactdefaultargs") NotifySubMeshesHypothesisModification;
		%feature("autodoc", "	* /*! * \brief Say all submeshes that theChangedHyp has been modified */

	:param theChangedHyp:
	:type theChangedHyp: SMESH_Hypothesis *
	:rtype: None
") NotifySubMeshesHypothesisModification;
		void NotifySubMeshesHypothesisModification (const SMESH_Hypothesis * theChangedHyp);
		%feature("compactdefaultargs") GetSubMeshUsingHypothesis;
		%feature("autodoc", "	:param anHyp:
	:type anHyp: SMESHDS_Hypothesis *
	:rtype: std::list< SMESH_subMesh *>
") GetSubMeshUsingHypothesis;
		const std::list< SMESH_subMesh *> & GetSubMeshUsingHypothesis (SMESHDS_Hypothesis * anHyp);
		%feature("compactdefaultargs") IsUsedHypothesis;
		%feature("autodoc", "	* /*! * \brief Return True if anHyp is used to mesh aSubShape */

	:param anHyp:
	:type anHyp: SMESHDS_Hypothesis *
	:param aSubMesh:
	:type aSubMesh: SMESH_subMesh *
	:rtype: bool
") IsUsedHypothesis;
		bool IsUsedHypothesis (SMESHDS_Hypothesis * anHyp,const SMESH_subMesh * aSubMesh);
		%feature("compactdefaultargs") IsNotConformAllowed;
		%feature("autodoc", "	* /*! * \brief check if a hypothesis alowing notconform mesh is present */

	:rtype: bool
") IsNotConformAllowed;
		bool IsNotConformAllowed ();
		%feature("compactdefaultargs") IsMainShape;
		%feature("autodoc", "	:param theShape:
	:type theShape: TopoDS_Shape &
	:rtype: bool
") IsMainShape;
		bool IsMainShape (const TopoDS_Shape & theShape);
		%feature("compactdefaultargs") GetAncestors;
		%feature("autodoc", "	* /*! * \brief Return list of ancestors of theSubShape in the order * that lower dimention shapes come first */

	:param theSubShape:
	:type theSubShape: TopoDS_Shape &
	:rtype: TopTools_ListOfShape
") GetAncestors;
		const TopTools_ListOfShape & GetAncestors (const TopoDS_Shape & theSubShape);
		%feature("compactdefaultargs") SetAutoColor;
		%feature("autodoc", "	:param theAutoColor:
	:type theAutoColor: bool
	:rtype: None
") SetAutoColor;
		void SetAutoColor (bool theAutoColor);
		%feature("compactdefaultargs") GetAutoColor;
		%feature("autodoc", "	:rtype: bool
") GetAutoColor;
		bool GetAutoColor ();
		%feature("compactdefaultargs") GetAncestorMap;
		%feature("autodoc", "	:rtype: TopTools_IndexedDataMapOfShapeListOfShape
") GetAncestorMap;
		const TopTools_IndexedDataMapOfShapeListOfShape & GetAncestorMap ();
		%feature("compactdefaultargs") HasDuplicatedGroupNamesMED;
		%feature("autodoc", "	* /*! * \brief Check group names for duplications. * Consider maximum group name length stored in MED file */

	:rtype: bool
") HasDuplicatedGroupNamesMED;
		bool HasDuplicatedGroupNamesMED ();
		%feature("compactdefaultargs") ExportMED;
		%feature("autodoc", "	:param file:
	:type file: char *
	:param theMeshName: default value is NULL
	:type theMeshName: char *
	:param theAutoGroups: default value is true
	:type theAutoGroups: bool
	:param theVersion: default value is 0
	:type theVersion: int
	:rtype: None
") ExportMED;
		void ExportMED (const char * file,const char * theMeshName = NULL,bool theAutoGroups = true,int theVersion = 0);
		%feature("compactdefaultargs") ExportDAT;
		%feature("autodoc", "	:param file:
	:type file: char *
	:rtype: None
") ExportDAT;
		void ExportDAT (const char * file);
		%feature("compactdefaultargs") ExportUNV;
		%feature("autodoc", "	:param file:
	:type file: char *
	:rtype: None
") ExportUNV;
		void ExportUNV (const char * file);
		%feature("compactdefaultargs") ExportSTL;
		%feature("autodoc", "	:param file:
	:type file: char *
	:param isascii:
	:type isascii: bool
	:rtype: None
") ExportSTL;
		void ExportSTL (const char * file,const bool isascii);
		%feature("compactdefaultargs") NbNodes;
		%feature("autodoc", "	:rtype: int
") NbNodes;
		int NbNodes ();
		%feature("compactdefaultargs") Nb0DElements;
		%feature("autodoc", "	:rtype: int
") Nb0DElements;
		int Nb0DElements ();
		%feature("compactdefaultargs") NbEdges;
		%feature("autodoc", "	:param order: default value is ORDER_ANY
	:type order: SMDSAbs_ElementOrder
	:rtype: int
") NbEdges;
		int NbEdges (SMDSAbs_ElementOrder order = ORDER_ANY);
		%feature("compactdefaultargs") NbFaces;
		%feature("autodoc", "	:param order: default value is ORDER_ANY
	:type order: SMDSAbs_ElementOrder
	:rtype: int
") NbFaces;
		int NbFaces (SMDSAbs_ElementOrder order = ORDER_ANY);
		%feature("compactdefaultargs") NbTriangles;
		%feature("autodoc", "	:param order: default value is ORDER_ANY
	:type order: SMDSAbs_ElementOrder
	:rtype: int
") NbTriangles;
		int NbTriangles (SMDSAbs_ElementOrder order = ORDER_ANY);
		%feature("compactdefaultargs") NbQuadrangles;
		%feature("autodoc", "	:param order: default value is ORDER_ANY
	:type order: SMDSAbs_ElementOrder
	:rtype: int
") NbQuadrangles;
		int NbQuadrangles (SMDSAbs_ElementOrder order = ORDER_ANY);
		%feature("compactdefaultargs") NbPolygons;
		%feature("autodoc", "	:rtype: int
") NbPolygons;
		int NbPolygons ();
		%feature("compactdefaultargs") NbVolumes;
		%feature("autodoc", "	:param order: default value is ORDER_ANY
	:type order: SMDSAbs_ElementOrder
	:rtype: int
") NbVolumes;
		int NbVolumes (SMDSAbs_ElementOrder order = ORDER_ANY);
		%feature("compactdefaultargs") NbTetras;
		%feature("autodoc", "	:param order: default value is ORDER_ANY
	:type order: SMDSAbs_ElementOrder
	:rtype: int
") NbTetras;
		int NbTetras (SMDSAbs_ElementOrder order = ORDER_ANY);
		%feature("compactdefaultargs") NbHexas;
		%feature("autodoc", "	:param order: default value is ORDER_ANY
	:type order: SMDSAbs_ElementOrder
	:rtype: int
") NbHexas;
		int NbHexas (SMDSAbs_ElementOrder order = ORDER_ANY);
		%feature("compactdefaultargs") NbPyramids;
		%feature("autodoc", "	:param order: default value is ORDER_ANY
	:type order: SMDSAbs_ElementOrder
	:rtype: int
") NbPyramids;
		int NbPyramids (SMDSAbs_ElementOrder order = ORDER_ANY);
		%feature("compactdefaultargs") NbPrisms;
		%feature("autodoc", "	:param order: default value is ORDER_ANY
	:type order: SMDSAbs_ElementOrder
	:rtype: int
") NbPrisms;
		int NbPrisms (SMDSAbs_ElementOrder order = ORDER_ANY);
		%feature("compactdefaultargs") NbPolyhedrons;
		%feature("autodoc", "	:rtype: int
") NbPolyhedrons;
		int NbPolyhedrons ();
		%feature("compactdefaultargs") NbSubMesh;
		%feature("autodoc", "	:rtype: int
") NbSubMesh;
		int NbSubMesh ();
		%feature("compactdefaultargs") NbGroup;
		%feature("autodoc", "	:rtype: int
") NbGroup;
		int NbGroup ();
		%feature("compactdefaultargs") AddGroup;
		%feature("autodoc", "	:param theType:
	:type theType: SMDSAbs_ElementType
	:param theName:
	:type theName: char *
	:param theId:
	:type theId: int &
	:param theShape: default value is TopoDS_Shape()
	:type theShape: TopoDS_Shape &
	:rtype: SMESH_Group *
") AddGroup;
		SMESH_Group * AddGroup (const SMDSAbs_ElementType theType,const char * theName,Standard_Integer &OutValue,const TopoDS_Shape & theShape = TopoDS_Shape());
		%feature("compactdefaultargs") GetGroupIds;
		%feature("autodoc", "	:rtype: std::list<int>
") GetGroupIds;
		std::list<int> GetGroupIds ();
		%feature("compactdefaultargs") RemoveGroup;
		%feature("autodoc", "	:param theGroupID:
	:type theGroupID: int
	:rtype: None
") RemoveGroup;
		void RemoveGroup (const int theGroupID);
		%feature("compactdefaultargs") ConvertToStandalone;
		%feature("autodoc", "	:param theGroupID:
	:type theGroupID: int
	:rtype: SMESH_Group *
") ConvertToStandalone;
		SMESH_Group * ConvertToStandalone (int theGroupID);
		%feature("compactdefaultargs") GetElementType;
		%feature("autodoc", "	:param id:
	:type id: int
	:param iselem:
	:type iselem: bool
	:rtype: SMDSAbs_ElementType
") GetElementType;
		SMDSAbs_ElementType GetElementType (const int id,const bool iselem);
		%feature("compactdefaultargs") Dump;
		%feature("autodoc", "	:param save:
	:type save: ostream &
	:rtype: ostream
") Dump;
		ostream & Dump (ostream & save);
};


%nodefaultctor SMESH_MeshEditor_PathPoint;
class SMESH_MeshEditor_PathPoint {
	public:
		%feature("compactdefaultargs") SMESH_MeshEditor_PathPoint;
		%feature("autodoc", "	:rtype: None
") SMESH_MeshEditor_PathPoint;
		 SMESH_MeshEditor_PathPoint ();
		%feature("compactdefaultargs") SetPnt;
		%feature("autodoc", "	:param aP3D:
	:type aP3D: gp_Pnt
	:rtype: None
") SetPnt;
		void SetPnt (const gp_Pnt & aP3D);
		%feature("compactdefaultargs") SetTangent;
		%feature("autodoc", "	:param aTgt:
	:type aTgt: gp_Dir
	:rtype: None
") SetTangent;
		void SetTangent (const gp_Dir & aTgt);
		%feature("compactdefaultargs") SetAngle;
		%feature("autodoc", "	:param aBeta:
	:type aBeta: double &
	:rtype: None
") SetAngle;
		void SetAngle (const double & aBeta);
		%feature("compactdefaultargs") SetParameter;
		%feature("autodoc", "	:param aPrm:
	:type aPrm: double &
	:rtype: None
") SetParameter;
		void SetParameter (const double & aPrm);
		%feature("compactdefaultargs") Pnt;
		%feature("autodoc", "	:rtype: gp_Pnt
") Pnt;
		const gp_Pnt  Pnt ();
		%feature("compactdefaultargs") Tangent;
		%feature("autodoc", "	:rtype: gp_Dir
") Tangent;
		const gp_Dir  Tangent ();
		%feature("compactdefaultargs") Angle;
		%feature("autodoc", "	:rtype: double
") Angle;
		double Angle ();
		%feature("compactdefaultargs") Parameter;
		%feature("autodoc", "	:rtype: double
") Parameter;
		double Parameter ();
};


%nodefaultctor SMESH_MeshVSLink;
class SMESH_MeshVSLink : public MeshVS_DataSource3D {
	public:
		%feature("compactdefaultargs") SMESH_MeshVSLink;
		%feature("autodoc", "	* Constructor

	:param aMesh:
	:type aMesh: SMESH_Mesh *
	:rtype: None
") SMESH_MeshVSLink;
		 SMESH_MeshVSLink (const SMESH_Mesh * aMesh);
		%feature("compactdefaultargs") GetGeom;
		%feature("autodoc", "	* Returns geometry information about node ( if IsElement is False ) or element ( IsElement is True ) by co-ordinates. For element this method must return all its nodes co-ordinates in the strict order: X, Y, Z and with nodes order is the same as in wire bounding the face or link. NbNodes is number of nodes of element. It is recommended to return 1 for node. Type is an element type.

	:param ID:
	:type ID: int
	:param IsElement:
	:type IsElement: bool
	:param Coords:
	:type Coords: TColStd_Array1OfReal &
	:param NbNodes:
	:type NbNodes: int &
	:param Type:
	:type Type: MeshVS_EntityType &
	:rtype: bool
") GetGeom;
		Standard_Boolean GetGeom (const Standard_Integer ID,const Standard_Boolean IsElement,TColStd_Array1OfReal & Coords,Standard_Integer &OutValue,MeshVS_EntityType & Type);
		%feature("compactdefaultargs") Get3DGeom;
		%feature("autodoc", "	:param ID:
	:type ID: int
	:param NbNodes:
	:type NbNodes: int &
	:param Data:
	:type Data: Handle_MeshVS_HArray1OfSequenceOfInteger &
	:rtype: bool
") Get3DGeom;
		Standard_Boolean Get3DGeom (const Standard_Integer ID,Standard_Integer &OutValue,Handle_MeshVS_HArray1OfSequenceOfInteger & Data);
		%feature("compactdefaultargs") GetGeomType;
		%feature("autodoc", "	* This method is similar to GetGeom, but returns only element or node type. This method is provided for a fine performance.

	:param ID:
	:type ID: int
	:param IsElement:
	:type IsElement: bool
	:param Type:
	:type Type: MeshVS_EntityType &
	:rtype: bool
") GetGeomType;
		Standard_Boolean GetGeomType (const Standard_Integer ID,const Standard_Boolean IsElement,MeshVS_EntityType & Type);
		%feature("compactdefaultargs") GetAddr;
		%feature("autodoc", "	* This method returns by number an address of any entity which represents element or node data structure.

	:param ID:
	:type ID: int
	:param IsElement:
	:type IsElement: bool
	:rtype: Standard_Address
") GetAddr;
		Standard_Address GetAddr (const Standard_Integer ID,const Standard_Boolean IsElement);
		%feature("compactdefaultargs") GetNodesByElement;
		%feature("autodoc", "	* This method returns information about what node this element consist of.

	:param ID:
	:type ID: int
	:param NodeIDs:
	:type NodeIDs: TColStd_Array1OfInteger &
	:param NbNodes:
	:type NbNodes: int &
	:rtype: bool
") GetNodesByElement;
		Standard_Boolean GetNodesByElement (const Standard_Integer ID,TColStd_Array1OfInteger & NodeIDs,Standard_Integer &OutValue);
		%feature("compactdefaultargs") GetAllNodes;
		%feature("autodoc", "	* This method returns map of all nodes the object consist of.

	:rtype: TColStd_PackedMapOfInteger
") GetAllNodes;
		const TColStd_PackedMapOfInteger & GetAllNodes ();
		%feature("compactdefaultargs") GetAllElements;
		%feature("autodoc", "	* This method returns map of all elements the object consist of.

	:rtype: TColStd_PackedMapOfInteger
") GetAllElements;
		const TColStd_PackedMapOfInteger & GetAllElements ();
		%feature("compactdefaultargs") GetNormal;
		%feature("autodoc", "	* This method calculates normal of face, which is using for correct reflection presentation. There is default method, for advance reflection this method can be redefined.

	:param Id:
	:type Id: int
	:param Max:
	:type Max: int
	:param nx:
	:type nx: float &
	:param ny:
	:type ny: float &
	:param nz:
	:type nz: float &
	:rtype: bool
") GetNormal;
		Standard_Boolean GetNormal (const Standard_Integer Id,const Standard_Integer Max,Standard_Real &OutValue,Standard_Real &OutValue,Standard_Real &OutValue);
		%feature("compactdefaultargs") GetAllGroups;
		%feature("autodoc", "	* This method returns map of all groups the object contains.

	:param Ids:
	:type Ids: TColStd_PackedMapOfInteger &
	:rtype: None
") GetAllGroups;
		void GetAllGroups (TColStd_PackedMapOfInteger & Ids);
		%feature("compactdefaultargs") DynamicType;
		%feature("autodoc", "	:rtype: Handle_Standard_Type
") DynamicType;
		Handle_Standard_Type DynamicType ();
};


%nodefaultctor SMESH_MesherHelper;
class SMESH_MesherHelper {
/* public enums */
enum MType {
	LINEAR = 0,
	QUADRATIC = 1,
	COMP = 2,
};

/* end public enums declaration */

	public:
		%feature("compactdefaultargs") IsMedium;
		%feature("autodoc", "	* /*! * \brief Returns true if given node is medium * \param n - node to check * \param typeToCheck - type of elements containing the node to ask about node status * etval bool - check result */

	:param node:
	:type node: SMDS_MeshNode *
	:param typeToCheck: default value is SMDSAbs_All
	:type typeToCheck: SMDSAbs_ElementType
	:rtype: bool
") IsMedium;
		static bool IsMedium (const SMDS_MeshNode * node,const SMDSAbs_ElementType typeToCheck = SMDSAbs_All);
		%feature("compactdefaultargs") LoadNodeColumns;
		%feature("autodoc", "	* /*! * \brief Load nodes bound to face into a map of node columns * \param theParam2ColumnMap - map of node columns to fill * \param theFace - the face on which nodes are searched for * \param theBaseEdge - the edge nodes of which are columns' bases * \param theMesh - the mesh containing nodes * etval bool - false if something is wrong * * The key of the map is a normalized parameter of each * base node on theBaseEdge. * This method works in supposition that nodes on the face * forms a rectangular grid and elements can be quardrangles or triangles */

	:param theParam2ColumnMap:
	:type theParam2ColumnMap: TParam2ColumnMap &
	:param theFace:
	:type theFace: TopoDS_Face &
	:param theBaseEdge:
	:type theBaseEdge: TopoDS_Edge &
	:param theMesh:
	:type theMesh: SMESHDS_Mesh *
	:rtype: bool
") LoadNodeColumns;
		static bool LoadNodeColumns (TParam2ColumnMap & theParam2ColumnMap,const TopoDS_Face & theFace,const TopoDS_Edge & theBaseEdge,SMESHDS_Mesh * theMesh);
		%feature("compactdefaultargs") GetSubShapeByNode;
		%feature("autodoc", "	* /*! * \brief Return support shape of a node * \param node - the node * \param meshDS - mesh DS * etval TopoDS_Shape - found support shape */

	:param node:
	:type node: SMDS_MeshNode *
	:param meshDS:
	:type meshDS: SMESHDS_Mesh *
	:rtype: TopoDS_Shape
") GetSubShapeByNode;
		static TopoDS_Shape GetSubShapeByNode (const SMDS_MeshNode * node,SMESHDS_Mesh * meshDS);
		%feature("compactdefaultargs") WrapIndex;
		%feature("autodoc", "	* /*! * \brief Return a valid node index, fixing the given one if necessary * \param ind - node index * \param nbNodes - total nb of nodes * etval int - valid node index */

	:param ind:
	:type ind: int
	:param nbNodes:
	:type nbNodes: int
	:rtype: int
") WrapIndex;
		static int WrapIndex (const int ind,const int nbNodes);
		%feature("compactdefaultargs") NbAncestors;
		%feature("autodoc", "	* /*! * \brief Return number of unique ancestors of the shape */

	:param shape:
	:type shape: TopoDS_Shape &
	:param mesh:
	:type mesh: SMESH_Mesh &
	:param ancestorType: default value is TopAbs_SHAPE
	:type ancestorType: TopAbs_ShapeEnum
	:rtype: int
") NbAncestors;
		static int NbAncestors (const TopoDS_Shape & shape,const SMESH_Mesh & mesh,TopAbs_ShapeEnum ancestorType = TopAbs_SHAPE);
		%feature("compactdefaultargs") SMESH_MesherHelper;
		%feature("autodoc", "	:param theMesh:
	:type theMesh: SMESH_Mesh &
	:rtype: None
") SMESH_MesherHelper;
		 SMESH_MesherHelper (SMESH_Mesh & theMesh);
		%feature("compactdefaultargs") GetMesh;
		%feature("autodoc", "	:rtype: SMESH_Mesh *
") GetMesh;
		SMESH_Mesh * GetMesh ();
		%feature("compactdefaultargs") GetMeshDS;
		%feature("autodoc", "	:rtype: SMESHDS_Mesh *
") GetMeshDS;
		SMESHDS_Mesh * GetMeshDS ();
		%feature("compactdefaultargs") IsQuadraticSubMesh;
		%feature("autodoc", "	* /*! * Check submesh for given shape: if all elements on this shape are quadratic, * quadratic elements will be created. Also fill myTLinkNodeMap */

	:param theShape:
	:type theShape: TopoDS_Shape &
	:rtype: bool
") IsQuadraticSubMesh;
		bool IsQuadraticSubMesh (const TopoDS_Shape & theShape);
		%feature("compactdefaultargs") SetIsQuadratic;
		%feature("autodoc", "	* /*! * \brief Set order of elements to create without calling IsQuadraticSubMesh() */

	:param theBuildQuadratic:
	:type theBuildQuadratic: bool
	:rtype: None
") SetIsQuadratic;
		void SetIsQuadratic (const bool theBuildQuadratic);
		%feature("compactdefaultargs") GetIsQuadratic;
		%feature("autodoc", "	* /*! * \brief Return myCreateQuadratic flag */

	:rtype: bool
") GetIsQuadratic;
		bool GetIsQuadratic ();
		%feature("compactdefaultargs") FixQuadraticElements;
		%feature("autodoc", "	* /*! * \brief Move medium nodes of faces and volumes to fix distorted elements * \param volumeOnly - fix nodes on geom faces or not if the shape is solid */

	:param volumeOnly: default value is true
	:type volumeOnly: bool
	:rtype: None
") FixQuadraticElements;
		void FixQuadraticElements (bool volumeOnly = true);
		%feature("compactdefaultargs") SetElementsOnShape;
		%feature("autodoc", "	* /*! * \brief To set created elements on the shape set by IsQuadraticSubMesh() * or the next methods. By defaul elements are set on the shape if * a mesh has no shape to be meshed */

	:param toSet:
	:type toSet: bool
	:rtype: None
") SetElementsOnShape;
		void SetElementsOnShape (bool toSet);
		%feature("compactdefaultargs") SetSubShape;
		%feature("autodoc", "	* /*! * \brief Set shape to make elements on without calling IsQuadraticSubMesh() */

	:param subShapeID:
	:type subShapeID: int
	:rtype: None
") SetSubShape;
		void SetSubShape (const int subShapeID);
		%feature("compactdefaultargs") SetSubShape;
		%feature("autodoc", "	* //!==SMESHDS_Mesh::ShapeToIndex(shape)

	:param subShape:
	:type subShape: TopoDS_Shape &
	:rtype: None
") SetSubShape;
		void SetSubShape (const TopoDS_Shape & subShape);
		%feature("compactdefaultargs") GetSubShapeID;
		%feature("autodoc", "	* /*! * \brief Return ID of the shape set by IsQuadraticSubMesh() or SetSubShape() * etval int - shape index in SMESHDS */

	:rtype: int
") GetSubShapeID;
		int GetSubShapeID ();
		%feature("compactdefaultargs") GetSubShape;
		%feature("autodoc", "	* /*! * \brief Return the shape set by IsQuadraticSubMesh() or SetSubShape() */

	:rtype: TopoDS_Shape
") GetSubShape;
		TopoDS_Shape GetSubShape ();
		%feature("compactdefaultargs") AddNode;
		%feature("autodoc", "	* /*! * Creates a node */

	:param x:
	:type x: double
	:param y:
	:type y: double
	:param z:
	:type z: double
	:param ID: default value is 0
	:type ID: int
	:rtype: SMDS_MeshNode *
") AddNode;
		SMDS_MeshNode * AddNode (double x,double y,double z,int ID = 0);
		%feature("compactdefaultargs") AddEdge;
		%feature("autodoc", "	* /*! * Creates quadratic or linear edge */

	:param n1:
	:type n1: SMDS_MeshNode *
	:param n2:
	:type n2: SMDS_MeshNode *
	:param id: default value is 0
	:type id: int
	:param force3d: default value is true
	:type force3d: bool
	:rtype: SMDS_MeshEdge *
") AddEdge;
		SMDS_MeshEdge * AddEdge (const SMDS_MeshNode * n1,const SMDS_MeshNode * n2,const int id = 0,const bool force3d = true);
		%feature("compactdefaultargs") AddFace;
		%feature("autodoc", "	* /*! * Creates quadratic or linear triangle */

	:param n1:
	:type n1: SMDS_MeshNode *
	:param n2:
	:type n2: SMDS_MeshNode *
	:param n3:
	:type n3: SMDS_MeshNode *
	:param id: default value is 0
	:type id: int
	:param force3d: default value is false
	:type force3d: bool
	:rtype: SMDS_MeshFace *
") AddFace;
		SMDS_MeshFace * AddFace (const SMDS_MeshNode * n1,const SMDS_MeshNode * n2,const SMDS_MeshNode * n3,const int id = 0,const bool force3d = false);
		%feature("compactdefaultargs") AddFace;
		%feature("autodoc", "	* /*! * Creates quadratic or linear quadrangle */

	:param n1:
	:type n1: SMDS_MeshNode *
	:param n2:
	:type n2: SMDS_MeshNode *
	:param n3:
	:type n3: SMDS_MeshNode *
	:param n4:
	:type n4: SMDS_MeshNode *
	:param id: default value is 0
	:type id: int
	:param force3d: default value is false
	:type force3d: bool
	:rtype: SMDS_MeshFace *
") AddFace;
		SMDS_MeshFace * AddFace (const SMDS_MeshNode * n1,const SMDS_MeshNode * n2,const SMDS_MeshNode * n3,const SMDS_MeshNode * n4,const int id = 0,const bool force3d = false);
		%feature("compactdefaultargs") AddVolume;
		%feature("autodoc", "	* /*! * Creates quadratic or linear tetraahedron */

	:param n1:
	:type n1: SMDS_MeshNode *
	:param n2:
	:type n2: SMDS_MeshNode *
	:param n3:
	:type n3: SMDS_MeshNode *
	:param n4:
	:type n4: SMDS_MeshNode *
	:param id: default value is 0
	:type id: int
	:param force3d: default value is true
	:type force3d: bool
	:rtype: SMDS_MeshVolume *
") AddVolume;
		SMDS_MeshVolume * AddVolume (const SMDS_MeshNode * n1,const SMDS_MeshNode * n2,const SMDS_MeshNode * n3,const SMDS_MeshNode * n4,const int id = 0,const bool force3d = true);
		%feature("compactdefaultargs") AddVolume;
		%feature("autodoc", "	* /*! * Creates quadratic or linear pyramid */

	:param n1:
	:type n1: SMDS_MeshNode *
	:param n2:
	:type n2: SMDS_MeshNode *
	:param n3:
	:type n3: SMDS_MeshNode *
	:param n4:
	:type n4: SMDS_MeshNode *
	:param n5:
	:type n5: SMDS_MeshNode *
	:param id: default value is 0
	:type id: int
	:param force3d: default value is true
	:type force3d: bool
	:rtype: SMDS_MeshVolume *
") AddVolume;
		SMDS_MeshVolume * AddVolume (const SMDS_MeshNode * n1,const SMDS_MeshNode * n2,const SMDS_MeshNode * n3,const SMDS_MeshNode * n4,const SMDS_MeshNode * n5,const int id = 0,const bool force3d = true);
		%feature("compactdefaultargs") AddVolume;
		%feature("autodoc", "	* /*! * Creates quadratic or linear pentahedron */

	:param n1:
	:type n1: SMDS_MeshNode *
	:param n2:
	:type n2: SMDS_MeshNode *
	:param n3:
	:type n3: SMDS_MeshNode *
	:param n4:
	:type n4: SMDS_MeshNode *
	:param n5:
	:type n5: SMDS_MeshNode *
	:param n6:
	:type n6: SMDS_MeshNode *
	:param id: default value is 0
	:type id: int
	:param force3d: default value is true
	:type force3d: bool
	:rtype: SMDS_MeshVolume *
") AddVolume;
		SMDS_MeshVolume * AddVolume (const SMDS_MeshNode * n1,const SMDS_MeshNode * n2,const SMDS_MeshNode * n3,const SMDS_MeshNode * n4,const SMDS_MeshNode * n5,const SMDS_MeshNode * n6,const int id = 0,const bool force3d = true);
		%feature("compactdefaultargs") AddVolume;
		%feature("autodoc", "	* /*! * Creates quadratic or linear hexahedron */

	:param n1:
	:type n1: SMDS_MeshNode *
	:param n2:
	:type n2: SMDS_MeshNode *
	:param n3:
	:type n3: SMDS_MeshNode *
	:param n4:
	:type n4: SMDS_MeshNode *
	:param n5:
	:type n5: SMDS_MeshNode *
	:param n6:
	:type n6: SMDS_MeshNode *
	:param n7:
	:type n7: SMDS_MeshNode *
	:param n8:
	:type n8: SMDS_MeshNode *
	:param id: default value is 0
	:type id: int
	:param force3d: default value is true
	:type force3d: bool
	:rtype: SMDS_MeshVolume *
") AddVolume;
		SMDS_MeshVolume * AddVolume (const SMDS_MeshNode * n1,const SMDS_MeshNode * n2,const SMDS_MeshNode * n3,const SMDS_MeshNode * n4,const SMDS_MeshNode * n5,const SMDS_MeshNode * n6,const SMDS_MeshNode * n7,const SMDS_MeshNode * n8,const int id = 0,bool force3d = true);
		%feature("compactdefaultargs") GetNodeU;
		%feature("autodoc", "	* /*! * \brief Return U of the given node on the edge */

	:param theEdge:
	:type theEdge: TopoDS_Edge &
	:param theNode:
	:type theNode: SMDS_MeshNode *
	:param check: default value is 0
	:type check: bool *
	:rtype: double
") GetNodeU;
		double GetNodeU (const TopoDS_Edge & theEdge,const SMDS_MeshNode * theNode,bool * check = 0);
		%feature("compactdefaultargs") GetNodeUV;
		%feature("autodoc", "	* /*! * \brief Return node UV on face * \param inFaceNode - a node of element being created located inside a face */

	:param F:
	:type F: TopoDS_Face &
	:param n:
	:type n: SMDS_MeshNode *
	:param inFaceNode: default value is 0
	:type inFaceNode: SMDS_MeshNode *
	:param check: default value is 0
	:type check: bool *
	:rtype: gp_XY
") GetNodeUV;
		gp_XY GetNodeUV (const TopoDS_Face & F,const SMDS_MeshNode * n,const SMDS_MeshNode * inFaceNode = 0,bool * check = 0);
		%feature("compactdefaultargs") CheckNodeUV;
		%feature("autodoc", "	* /*! * \brief Check and fix node UV on a face * etval bool - false if UV is bad and could not be fixed */

	:param F:
	:type F: TopoDS_Face &
	:param n:
	:type n: SMDS_MeshNode *
	:param uv:
	:type uv: gp_XY
	:param tol:
	:type tol: double
	:rtype: bool
") CheckNodeUV;
		bool CheckNodeUV (const TopoDS_Face & F,const SMDS_MeshNode * n,gp_XY & uv,const double tol);
		%feature("compactdefaultargs") GetMiddleUV;
		%feature("autodoc", "	* /*! * \brief Return middle UV taking in account surface period */

	:param surface:
	:type surface: Handle_Geom_Surface &
	:param uv1:
	:type uv1: gp_XY
	:param uv2:
	:type uv2: gp_XY
	:rtype: gp_XY
") GetMiddleUV;
		static gp_XY GetMiddleUV (const Handle_Geom_Surface & surface,const gp_XY & uv1,const gp_XY & uv2);
		%feature("compactdefaultargs") GetNodeUVneedInFaceNode;
		%feature("autodoc", "	* /*! * \brief Check if inFaceNode argument is necessary for call GetNodeUV(F,..) * etval bool - return true if the face is periodic * * if F is Null, answer about subshape set through IsQuadraticSubMesh() or * SetSubShape() */

	:param F: default value is TopoDS_Face()
	:type F: TopoDS_Face &
	:rtype: bool
") GetNodeUVneedInFaceNode;
		bool GetNodeUVneedInFaceNode (const TopoDS_Face & F = TopoDS_Face());
		%feature("compactdefaultargs") IsDegenShape;
		%feature("autodoc", "	* /*! * \brief Check if shape is a degenerated edge or it's vertex * \param subShape - edge or vertex index in SMESHDS * etval bool - true if subShape is a degenerated shape * * It works only if IsQuadraticSubMesh() or SetSubShape() has been called */

	:param subShape:
	:type subShape: int
	:rtype: bool
") IsDegenShape;
		bool IsDegenShape (const int subShape);
		%feature("compactdefaultargs") IsSeamShape;
		%feature("autodoc", "	* /*! * \brief Check if shape is a seam edge or it's vertex * \param subShape - edge or vertex index in SMESHDS * etval bool - true if subShape is a seam shape * * It works only if IsQuadraticSubMesh() or SetSubShape() has been called. * Seam shape has two 2D alternative represenations on the face */

	:param subShape:
	:type subShape: int
	:rtype: bool
") IsSeamShape;
		bool IsSeamShape (const int subShape);
		%feature("compactdefaultargs") IsSeamShape;
		%feature("autodoc", "	* /*! * \brief Check if shape is a seam edge or it's vertex * \param subShape - edge or vertex * etval bool - true if subShape is a seam shape * * It works only if IsQuadraticSubMesh() or SetSubShape() has been called. * Seam shape has two 2D alternative represenations on the face */

	:param subShape:
	:type subShape: TopoDS_Shape &
	:rtype: bool
") IsSeamShape;
		bool IsSeamShape (const TopoDS_Shape & subShape);
		%feature("compactdefaultargs") IsRealSeam;
		%feature("autodoc", "	* /*! * \brief Return true if an edge or a vertex encounters twice in face wire * \param subShape - Id of edge or vertex */

	:param subShape:
	:type subShape: int
	:rtype: bool
") IsRealSeam;
		bool IsRealSeam (const int subShape);
		%feature("compactdefaultargs") IsRealSeam;
		%feature("autodoc", "	* /*! * \brief Return true if an edge or a vertex encounters twice in face wire * \param subShape - edge or vertex */

	:param subShape:
	:type subShape: TopoDS_Shape &
	:rtype: bool
") IsRealSeam;
		bool IsRealSeam (const TopoDS_Shape & subShape);
		%feature("compactdefaultargs") HasSeam;
		%feature("autodoc", "	* /*! * \brief Check if the shape set through IsQuadraticSubMesh() or SetSubShape() * has a seam edge * etval bool - true if it has */

	:rtype: bool
") HasSeam;
		bool HasSeam ();
		%feature("compactdefaultargs") GetPeriodicIndex;
		%feature("autodoc", "	* /*! * \brief Return index of periodic parametric direction of a closed face * etval int - 1 for U, 2 for V direction */

	:rtype: int
") GetPeriodicIndex;
		int GetPeriodicIndex ();
		%feature("compactdefaultargs") GetOtherParam;
		%feature("autodoc", "	* /*! * \brief Return an alternative parameter for a node on seam */

	:param param:
	:type param: double
	:rtype: double
") GetOtherParam;
		double GetOtherParam (const double param);
		%feature("compactdefaultargs") GetMediumNode;
		%feature("autodoc", "	* /** * Special function for search or creation medium node */

	:param n1:
	:type n1: SMDS_MeshNode *
	:param n2:
	:type n2: SMDS_MeshNode *
	:param force3d:
	:type force3d: bool
	:rtype: SMDS_MeshNode *
") GetMediumNode;
		const SMDS_MeshNode * GetMediumNode (const SMDS_MeshNode * n1,const SMDS_MeshNode * n2,const bool force3d);
		%feature("compactdefaultargs") AddTLinkNode;
		%feature("autodoc", "	* /*! * Auxilary function for filling myTLinkNodeMap */

	:param n1:
	:type n1: SMDS_MeshNode *
	:param n2:
	:type n2: SMDS_MeshNode *
	:param n12:
	:type n12: SMDS_MeshNode *
	:rtype: None
") AddTLinkNode;
		void AddTLinkNode (const SMDS_MeshNode * n1,const SMDS_MeshNode * n2,const SMDS_MeshNode * n12);
		%feature("compactdefaultargs") AddTLinkNodeMap;
		%feature("autodoc", "	* /** * Auxilary function for filling myTLinkNodeMap */

	:param aMap:
	:type aMap: TLinkNodeMap &
	:rtype: None
") AddTLinkNodeMap;
		void AddTLinkNodeMap (const TLinkNodeMap & aMap);
		%feature("compactdefaultargs") GetTLinkNodeMap;
		%feature("autodoc", "	* /** * Returns myTLinkNodeMap */

	:rtype: TLinkNodeMap
") GetTLinkNodeMap;
		const TLinkNodeMap & GetTLinkNodeMap ();
};


%nodefaultctor SMESH_NodeSearcher;
class SMESH_NodeSearcher {
	public:
		%feature("compactdefaultargs") FindClosestTo;
		%feature("autodoc", "	:param pnt:
	:type pnt: gp_Pnt
	:rtype: SMDS_MeshNode *
") FindClosestTo;
		const SMDS_MeshNode * FindClosestTo (const gp_Pnt & pnt);
		%feature("compactdefaultargs") MoveNode;
		%feature("autodoc", "	:param node:
	:type node: SMDS_MeshNode *
	:param toPnt:
	:type toPnt: gp_Pnt
	:rtype: None
") MoveNode;
		void MoveNode (const SMDS_MeshNode * node,const gp_Pnt & toPnt);
};


%nodefaultctor SMESH_Pattern;
class SMESH_Pattern {
/* public enums */
enum ErrorCode {
	ERR_OK = 0,
	ERR_READ_NB_POINTS = 1,
	ERR_READ_POINT_COORDS = 2,
	ERR_READ_TOO_FEW_POINTS = 3,
	ERR_READ_3D_COORD = 4,
	ERR_READ_NO_KEYPOINT = 5,
	ERR_READ_BAD_INDEX = 6,
	ERR_READ_ELEM_POINTS = 7,
	ERR_READ_NO_ELEMS = 8,
	ERR_READ_BAD_KEY_POINT = 9,
	ERR_SAVE_NOT_LOADED = 10,
	ERR_LOAD_EMPTY_SUBMESH = 11,
	ERR_LOADF_NARROW_FACE = 12,
	ERR_LOADF_CLOSED_FACE = 13,
	ERR_LOADF_CANT_PROJECT = 14,
	ERR_LOADV_BAD_SHAPE = 15,
	ERR_LOADV_COMPUTE_PARAMS = 16,
	ERR_APPL_NOT_COMPUTED = 17,
	ERR_APPL_NOT_LOADED = 18,
	ERR_APPL_BAD_DIMENTION = 19,
	ERR_APPL_BAD_NB_VERTICES = 20,
	ERR_APPLF_BAD_TOPOLOGY = 21,
	ERR_APPLF_BAD_VERTEX = 22,
	ERR_APPLF_INTERNAL_EEROR = 23,
	ERR_APPLV_BAD_SHAPE = 24,
	ERR_APPLF_BAD_FACE_GEOM = 25,
	ERR_MAKEM_NOT_COMPUTED = 26,
};

/* end public enums declaration */

	public:
		%feature("compactdefaultargs") SMESH_Pattern;
		%feature("autodoc", "	:rtype: None
") SMESH_Pattern;
		 SMESH_Pattern ();
		%feature("compactdefaultargs") Clear;
		%feature("autodoc", "	:rtype: None
") Clear;
		void Clear ();
		%feature("compactdefaultargs") Load;
		%feature("autodoc", "	:param theFileContents:
	:type theFileContents: char *
	:rtype: bool
") Load;
		bool Load (const char * theFileContents);
		%feature("compactdefaultargs") Load;
		%feature("autodoc", "	:param theMesh:
	:type theMesh: SMESH_Mesh *
	:param theFace:
	:type theFace: TopoDS_Face &
	:param theProject: default value is false
	:type theProject: bool
	:rtype: bool
") Load;
		bool Load (SMESH_Mesh * theMesh,const TopoDS_Face & theFace,bool theProject = false);
		%feature("compactdefaultargs") Load;
		%feature("autodoc", "	:param theMesh:
	:type theMesh: SMESH_Mesh *
	:param theBlock:
	:type theBlock: TopoDS_Shell &
	:rtype: bool
") Load;
		bool Load (SMESH_Mesh * theMesh,const TopoDS_Shell & theBlock);
		%feature("compactdefaultargs") Save;
		%feature("autodoc", "	:param theFile:
	:type theFile: std::ostream &
	:rtype: bool
") Save;
		bool Save (std::ostream & theFile);
		%feature("compactdefaultargs") Apply;
		%feature("autodoc", "	:param theFace:
	:type theFace: TopoDS_Face &
	:param theVertexOnKeyPoint1:
	:type theVertexOnKeyPoint1: TopoDS_Vertex &
	:param theReverse:
	:type theReverse: bool
	:rtype: bool
") Apply;
		bool Apply (const TopoDS_Face & theFace,const TopoDS_Vertex & theVertexOnKeyPoint1,const bool theReverse);
		%feature("compactdefaultargs") Apply;
		%feature("autodoc", "	:param theBlock:
	:type theBlock: TopoDS_Shell &
	:param theVertex000:
	:type theVertex000: TopoDS_Vertex &
	:param theVertex001:
	:type theVertex001: TopoDS_Vertex &
	:rtype: bool
") Apply;
		bool Apply (const TopoDS_Shell & theBlock,const TopoDS_Vertex & theVertex000,const TopoDS_Vertex & theVertex001);
		%feature("compactdefaultargs") Apply;
		%feature("autodoc", "	:param theFace:
	:type theFace: SMDS_MeshFace *
	:param theNodeIndexOnKeyPoint1:
	:type theNodeIndexOnKeyPoint1: int
	:param theReverse:
	:type theReverse: bool
	:rtype: bool
") Apply;
		bool Apply (const SMDS_MeshFace * theFace,const int theNodeIndexOnKeyPoint1,const bool theReverse);
		%feature("compactdefaultargs") Apply;
		%feature("autodoc", "	:param theMesh:
	:type theMesh: SMESH_Mesh *
	:param theFace:
	:type theFace: SMDS_MeshFace *
	:param theSurface:
	:type theSurface: TopoDS_Shape &
	:param theNodeIndexOnKeyPoint1:
	:type theNodeIndexOnKeyPoint1: int
	:param theReverse:
	:type theReverse: bool
	:rtype: bool
") Apply;
		bool Apply (SMESH_Mesh * theMesh,const SMDS_MeshFace * theFace,const TopoDS_Shape & theSurface,const int theNodeIndexOnKeyPoint1,const bool theReverse);
		%feature("compactdefaultargs") Apply;
		%feature("autodoc", "	:param theMesh:
	:type theMesh: SMESH_Mesh *
	:param theFaces:
	:type theFaces: std::set< SMDS_MeshFace *> &
	:param theNodeIndexOnKeyPoint1:
	:type theNodeIndexOnKeyPoint1: int
	:param theReverse:
	:type theReverse: bool
	:rtype: bool
") Apply;
		bool Apply (SMESH_Mesh * theMesh,std::set<const SMDS_MeshFace *> & theFaces,const int theNodeIndexOnKeyPoint1,const bool theReverse);
		%feature("compactdefaultargs") Apply;
		%feature("autodoc", "	:param theVolume:
	:type theVolume: SMDS_MeshVolume *
	:param theNode000Index:
	:type theNode000Index: int
	:param theNode001Index:
	:type theNode001Index: int
	:rtype: bool
") Apply;
		bool Apply (const SMDS_MeshVolume * theVolume,const int theNode000Index,const int theNode001Index);
		%feature("compactdefaultargs") Apply;
		%feature("autodoc", "	:param theVolumes:
	:type theVolumes: std::set< SMDS_MeshVolume *> &
	:param theNode000Index:
	:type theNode000Index: int
	:param theNode001Index:
	:type theNode001Index: int
	:rtype: bool
") Apply;
		bool Apply (std::set<const SMDS_MeshVolume *> & theVolumes,const int theNode000Index,const int theNode001Index);
		%feature("compactdefaultargs") GetMappedPoints;
		%feature("autodoc", "	:param thePoints:
	:type thePoints: std::list< gp_XYZ *>
	:rtype: bool
") GetMappedPoints;
		bool GetMappedPoints (std::list<const gp_XYZ *> & thePoints);
		%feature("compactdefaultargs") MakeMesh;
		%feature("autodoc", "	:param theMesh:
	:type theMesh: SMESH_Mesh *
	:param toCreatePolygons: default value is false
	:type toCreatePolygons: bool
	:param toCreatePolyedrs: default value is false
	:type toCreatePolyedrs: bool
	:rtype: bool
") MakeMesh;
		bool MakeMesh (SMESH_Mesh * theMesh,const bool toCreatePolygons = false,const bool toCreatePolyedrs = false);
		%feature("compactdefaultargs") GetErrorCode;
		%feature("autodoc", "	:rtype: SMESH_Pattern::ErrorCode
") GetErrorCode;
		SMESH_Pattern::ErrorCode GetErrorCode ();
		%feature("compactdefaultargs") IsLoaded;
		%feature("autodoc", "	:rtype: bool
") IsLoaded;
		bool IsLoaded ();
		%feature("compactdefaultargs") Is2D;
		%feature("autodoc", "	:rtype: bool
") Is2D;
		bool Is2D ();
		%feature("compactdefaultargs") GetPoints;
		%feature("autodoc", "	:param thePoints:
	:type thePoints: std::list< gp_XYZ *>
	:rtype: bool
") GetPoints;
		bool GetPoints (std::list<const gp_XYZ *> & thePoints);
		%feature("compactdefaultargs") GetKeyPointIDs;
		%feature("autodoc", "	:rtype: std::list< int>
") GetKeyPointIDs;
		const std::list< int> & GetKeyPointIDs ();
		%feature("compactdefaultargs") GetElementPointIDs;
		%feature("autodoc", "	:param applied:
	:type applied: bool
	:rtype: std::list< std::list< int> >
") GetElementPointIDs;
		const std::list< std::list< int> > & GetElementPointIDs (bool applied);
		%feature("compactdefaultargs") DumpPoints;
		%feature("autodoc", "	:rtype: None
") DumpPoints;
		void DumpPoints ();
		%feature("compactdefaultargs") GetSubShape;
		%feature("autodoc", "	:param i:
	:type i: int
	:rtype: TopoDS_Shape
") GetSubShape;
		TopoDS_Shape GetSubShape (const int i);
};


%nodefaultctor SMESH_TLink;
class SMESH_TLink : public NLink {
	public:
		%feature("compactdefaultargs") SMESH_TLink;
		%feature("autodoc", "	:param n1:
	:type n1: SMDS_MeshNode *
	:param n2:
	:type n2: SMDS_MeshNode *
	:rtype: None
") SMESH_TLink;
		 SMESH_TLink (const SMDS_MeshNode * n1,const SMDS_MeshNode * n2);
		%feature("compactdefaultargs") SMESH_TLink;
		%feature("autodoc", "	:param link:
	:type link: NLink &
	:rtype: None
") SMESH_TLink;
		 SMESH_TLink (const NLink & link);
		%feature("compactdefaultargs") node1;
		%feature("autodoc", "	:rtype: SMDS_MeshNode *
") node1;
		const SMDS_MeshNode * node1 ();
		%feature("compactdefaultargs") node2;
		%feature("autodoc", "	:rtype: SMDS_MeshNode *
") node2;
		const SMDS_MeshNode * node2 ();
};


%nodefaultctor SMESH_subMesh;
class SMESH_subMesh {
/* public enums */
enum compute_state {
	NOT_READY = 0,
	READY_TO_COMPUTE = 1,
	COMPUTE_OK = 2,
	FAILED_TO_COMPUTE = 3,
};

enum algo_state {
	NO_ALGO = 0,
	MISSING_HYP = 1,
	HYP_OK = 2,
};

enum algo_event {
	ADD_HYP = 0,
	ADD_ALGO = 1,
	REMOVE_HYP = 2,
	REMOVE_ALGO = 3,
	ADD_FATHER_HYP = 4,
	ADD_FATHER_ALGO = 5,
	REMOVE_FATHER_HYP = 6,
	REMOVE_FATHER_ALGO = 7,
	MODIF_HYP = 8,
};

enum compute_event {
	MODIF_ALGO_STATE = 0,
	COMPUTE = 1,
	CLEAN = 2,
	SUBMESH_COMPUTED = 3,
	SUBMESH_RESTORED = 4,
	MESH_ENTITY_REMOVED = 5,
	CHECK_COMPUTE_STATE = 6,
};

enum event_type {
	ALGO_EVENT = 0,
	COMPUTE_EVENT = 1,
};

/* end public enums declaration */

	public:
		%feature("compactdefaultargs") SMESH_subMesh;
		%feature("autodoc", "	:param Id:
	:type Id: int
	:param father:
	:type father: SMESH_Mesh *
	:param meshDS:
	:type meshDS: SMESHDS_Mesh *
	:param aSubShape:
	:type aSubShape: TopoDS_Shape &
	:rtype: None
") SMESH_subMesh;
		 SMESH_subMesh (int Id,SMESH_Mesh * father,SMESHDS_Mesh * meshDS,const TopoDS_Shape & aSubShape);
		%feature("compactdefaultargs") GetId;
		%feature("autodoc", "	:rtype: int
") GetId;
		int GetId ();
		%feature("compactdefaultargs") GetFather;
		%feature("autodoc", "	:rtype: SMESH_Mesh *
") GetFather;
		SMESH_Mesh * GetFather ();
		%feature("compactdefaultargs") GetSubMeshDS;
		%feature("autodoc", "	:rtype: SMESHDS_SubMesh *
") GetSubMeshDS;
		SMESHDS_SubMesh * GetSubMeshDS ();
		%feature("compactdefaultargs") CreateSubMeshDS;
		%feature("autodoc", "	:rtype: SMESHDS_SubMesh *
") CreateSubMeshDS;
		SMESHDS_SubMesh * CreateSubMeshDS ();
		%feature("compactdefaultargs") GetFirstToCompute;
		%feature("autodoc", "	:rtype: SMESH_subMesh *
") GetFirstToCompute;
		SMESH_subMesh * GetFirstToCompute ();
		%feature("compactdefaultargs") DependsOn;
		%feature("autodoc", "	:rtype: std::map< int, SMESH_subMesh *>
") DependsOn;
		const std::map< int, SMESH_subMesh *> & DependsOn ();
		%feature("compactdefaultargs") getDependsOnIterator;
		%feature("autodoc", "	* /*! * \brief Return iterator on the submeshes this one depends on */

	:param includeSelf:
	:type includeSelf: bool
	:param complexShapeFirst:
	:type complexShapeFirst: bool
	:rtype: SMESH_subMeshIteratorPtr
") getDependsOnIterator;
		SMESH_subMeshIteratorPtr getDependsOnIterator (const bool includeSelf,const bool complexShapeFirst);
		%feature("compactdefaultargs") GetSubShape;
		%feature("autodoc", "	:rtype: TopoDS_Shape
") GetSubShape;
		const TopoDS_Shape  GetSubShape ();
		%feature("compactdefaultargs") SetEventListener;
		%feature("autodoc", "	* /*! * \brief Sets an event listener and its data to a submesh * \param listener - the listener to store * \param data - the listener data to store * \param where - the submesh to store the listener and it's data * * The method remembers the submesh \awhere it puts the listener in order to delete * them when HYP_OK algo_state is lost * After being set, event listener is notified on each event of \awhere submesh. */

	:param listener:
	:type listener: EventListener *
	:param data:
	:type data: EventListenerData *
	:param where:
	:type where: SMESH_subMesh *
	:rtype: None
") SetEventListener;
		void SetEventListener (EventListener * listener,EventListenerData * data,SMESH_subMesh * where);
		%feature("compactdefaultargs") GetEventListenerData;
		%feature("autodoc", "	* /*! * \brief Return an event listener data * \param listener - the listener whose data is * etval EventListenerData* - found data, maybe NULL */

	:param listener:
	:type listener: EventListener *
	:rtype: EventListenerData *
") GetEventListenerData;
		EventListenerData * GetEventListenerData (EventListener * listener);
		%feature("compactdefaultargs") DeleteEventListener;
		%feature("autodoc", "	* /*! * \brief Unregister the listener and delete it and it's data * \param listener - the event listener to delete */

	:param listener:
	:type listener: EventListener *
	:rtype: None
") DeleteEventListener;
		void DeleteEventListener (EventListener * listener);
		%feature("compactdefaultargs") AlgoStateEngine;
		%feature("autodoc", "	:param event:
	:type event: int
	:param anHyp:
	:type anHyp: SMESH_Hypothesis *
	:rtype: SMESH_Hypothesis::Hypothesis_Status
") AlgoStateEngine;
		SMESH_Hypothesis::Hypothesis_Status AlgoStateEngine (int event,SMESH_Hypothesis * anHyp);
		%feature("compactdefaultargs") SubMeshesAlgoStateEngine;
		%feature("autodoc", "	:param event:
	:type event: int
	:param anHyp:
	:type anHyp: SMESH_Hypothesis *
	:rtype: SMESH_Hypothesis::Hypothesis_Status
") SubMeshesAlgoStateEngine;
		SMESH_Hypothesis::Hypothesis_Status SubMeshesAlgoStateEngine (int event,SMESH_Hypothesis * anHyp);
		%feature("compactdefaultargs") GetAlgoState;
		%feature("autodoc", "	:rtype: int
") GetAlgoState;
		int GetAlgoState ();
		%feature("compactdefaultargs") GetComputeState;
		%feature("autodoc", "	:rtype: int
") GetComputeState;
		int GetComputeState ();
		%feature("compactdefaultargs") GetComputeError;
		%feature("autodoc", "	:rtype: SMESH_ComputeErrorPtr
") GetComputeError;
		SMESH_ComputeErrorPtr & GetComputeError ();
		%feature("compactdefaultargs") DumpAlgoState;
		%feature("autodoc", "	:param isMain:
	:type isMain: bool
	:rtype: None
") DumpAlgoState;
		void DumpAlgoState (bool isMain);
		%feature("compactdefaultargs") ComputeStateEngine;
		%feature("autodoc", "	:param event:
	:type event: int
	:rtype: bool
") ComputeStateEngine;
		bool ComputeStateEngine (int event);
		%feature("compactdefaultargs") Evaluate;
		%feature("autodoc", "	:param aResMap:
	:type aResMap: MapShapeNbElems &
	:rtype: bool
") Evaluate;
		bool Evaluate (MapShapeNbElems & aResMap);
		%feature("compactdefaultargs") IsConform;
		%feature("autodoc", "	:param theAlgo:
	:type theAlgo: SMESH_Algo *
	:rtype: bool
") IsConform;
		bool IsConform (const SMESH_Algo * theAlgo);
		%feature("compactdefaultargs") CanAddHypothesis;
		%feature("autodoc", "	:param theHypothesis:
	:type theHypothesis: SMESH_Hypothesis *
	:rtype: bool
") CanAddHypothesis;
		bool CanAddHypothesis (const SMESH_Hypothesis * theHypothesis);
		%feature("compactdefaultargs") IsApplicableHypotesis;
		%feature("autodoc", "	:param theHypothesis:
	:type theHypothesis: SMESH_Hypothesis *
	:param theShapeType:
	:type theShapeType: TopAbs_ShapeEnum
	:rtype: bool
") IsApplicableHypotesis;
		static bool IsApplicableHypotesis (const SMESH_Hypothesis * theHypothesis,const TopAbs_ShapeEnum theShapeType);
		%feature("compactdefaultargs") IsApplicableHypotesis;
		%feature("autodoc", "	:param theHypothesis:
	:type theHypothesis: SMESH_Hypothesis *
	:rtype: bool
") IsApplicableHypotesis;
		bool IsApplicableHypotesis (const SMESH_Hypothesis * theHypothesis);
		%feature("compactdefaultargs") CheckConcurentHypothesis;
		%feature("autodoc", "	:param theHypType:
	:type theHypType: int
	:rtype: SMESH_Hypothesis::Hypothesis_Status
") CheckConcurentHypothesis;
		SMESH_Hypothesis::Hypothesis_Status CheckConcurentHypothesis (const int theHypType);
		%feature("compactdefaultargs") IsEmpty;
		%feature("autodoc", "	* /*! * \brief Return true if no mesh entities is bound to the submesh */

	:rtype: bool
") IsEmpty;
		bool IsEmpty ();
		%feature("compactdefaultargs") IsMeshComputed;
		%feature("autodoc", "	:rtype: bool
") IsMeshComputed;
		bool IsMeshComputed ();
		%feature("compactdefaultargs") SetIsAlwaysComputed;
		%feature("autodoc", "	* /*! * \brief Allow algo->Compute() if a subshape of lower dim is meshed but * none mesh entity is bound to it */

	:param isAlCo:
	:type isAlCo: bool
	:rtype: None
") SetIsAlwaysComputed;
		void SetIsAlwaysComputed (bool isAlCo);
		%feature("compactdefaultargs") IsAlwaysComputed;
		%feature("autodoc", "	:rtype: bool
") IsAlwaysComputed;
		bool IsAlwaysComputed ();
};


%nodefaultctor SMESH_subMeshEventListener;
class SMESH_subMeshEventListener {
	public:
		%feature("compactdefaultargs") SMESH_subMeshEventListener;
		%feature("autodoc", "	* //!< if true, it will be deleted by SMESH_subMesh

	:param isDeletable:
	:type isDeletable: bool
	:rtype: None
") SMESH_subMeshEventListener;
		 SMESH_subMeshEventListener (bool isDeletable);
		%feature("compactdefaultargs") IsDeletable;
		%feature("autodoc", "	:rtype: bool
") IsDeletable;
		bool IsDeletable ();
		%feature("compactdefaultargs") ProcessEvent;
		%feature("autodoc", "	* /*! * \brief Do something on a certain event * \param event - algo_event or compute_event itself (of SMESH_subMesh) * \param eventType - ALGO_EVENT or COMPUTE_EVENT (of SMESH_subMesh) * \param subMesh - the submesh where the event occures * \param data - listener data stored in the subMesh * \param hyp - hypothesis, if eventType is algo_event * * The base implementation translates CLEAN event to the subMesh stored * in the listener data. Also it sends SUBMESH_COMPUTED event in case of * successful COMPUTE event. */

	:param event:
	:type event: int
	:param eventType:
	:type eventType: int
	:param subMesh:
	:type subMesh: SMESH_subMesh *
	:param data:
	:type data: SMESH_subMeshEventListenerData *
	:param hyp: default value is 0
	:type hyp: SMESH_Hypothesis *
	:rtype: None
") ProcessEvent;
		void ProcessEvent (const int event,const int eventType,SMESH_subMesh * subMesh,SMESH_subMeshEventListenerData * data,const SMESH_Hypothesis * hyp = 0);
};


%nodefaultctor SMESH_subMeshEventListenerData;
class SMESH_subMeshEventListenerData {
	public:
		%feature("compactdefaultargs") SMESH_subMeshEventListenerData;
		%feature("autodoc", "	* //!< generally: submeshes depending

	:param isDeletable:
	:type isDeletable: bool
	:rtype: None
") SMESH_subMeshEventListenerData;
		 SMESH_subMeshEventListenerData (bool isDeletable);
		%feature("compactdefaultargs") IsDeletable;
		%feature("autodoc", "	:rtype: bool
") IsDeletable;
		bool IsDeletable ();
		%feature("compactdefaultargs") MakeData;
		%feature("autodoc", "	* /*! * \brief Create a default listener data. * \param dependentSM - subMesh to store * \param type - data type * etval SMESH_subMeshEventListenerData* - a new listener data * * See SMESH_subMeshEventListener::ProcessEvent() to know how the default * listener uses it (implementation is in SMESH_subMesh.cxx) */

	:param dependentSM:
	:type dependentSM: SMESH_subMesh *
	:param type: default value is 0
	:type type: int
	:rtype: SMESH_subMeshEventListenerData *
") MakeData;
		static SMESH_subMeshEventListenerData * MakeData (SMESH_subMesh * dependentSM,const int type = 0);
};


%nodefaultctor SMESH_Algo;
class SMESH_Algo : public SMESH_Hypothesis {
	public:
		%feature("compactdefaultargs") SaveTo;
		%feature("autodoc", "	* /*! * \brief Saves nothing in a stream * \param save - the stream * etval virtual std::ostream & - the stream */

	:param save:
	:type save: std::ostream &
	:rtype: std::ostream
") SaveTo;
		std::ostream & SaveTo (std::ostream & save);
		%feature("compactdefaultargs") LoadFrom;
		%feature("autodoc", "	* /*! * \brief Loads nothing from a stream * \param load - the stream * etval virtual std::ostream & - the stream */

	:param load:
	:type load: std::istream &
	:rtype: std::istream
") LoadFrom;
		std::istream & LoadFrom (std::istream & load);
		%feature("compactdefaultargs") GetCompatibleHypothesis;
		%feature("autodoc", "	* /*! * \brief Returns all types of compatible hypotheses */

	:rtype: std::vector< std::string>
") GetCompatibleHypothesis;
		const std::vector< std::string> & GetCompatibleHypothesis ();
		%feature("compactdefaultargs") CheckHypothesis;
		%feature("autodoc", "	* /*! * \brief Check hypothesis definition to mesh a shape * \param aMesh - the mesh * \param aShape - the shape * \param aStatus - check result * etval bool - true if hypothesis is well defined */

	:param aMesh:
	:type aMesh: SMESH_Mesh &
	:param aShape:
	:type aShape: TopoDS_Shape &
	:param aStatus:
	:type aStatus: SMESH_Hypothesis::Hypothesis_Status &
	:rtype: bool
") CheckHypothesis;
		bool CheckHypothesis (SMESH_Mesh & aMesh,const TopoDS_Shape & aShape,SMESH_Hypothesis::Hypothesis_Status & aStatus);
		%feature("compactdefaultargs") Compute;
		%feature("autodoc", "	* /*! * \brief Computes mesh on a shape * \param aMesh - the mesh * \param aShape - the shape * etval bool - is a success * * Algorithms that !NeedDescretBoundary() || !OnlyUnaryInput() are * to set SMESH_ComputeError returned by SMESH_submesh::GetComputeError() * to report problematic subshapes */

	:param aMesh:
	:type aMesh: SMESH_Mesh &
	:param aShape:
	:type aShape: TopoDS_Shape &
	:rtype: bool
") Compute;
		bool Compute (SMESH_Mesh & aMesh,const TopoDS_Shape & aShape);
		%feature("compactdefaultargs") Compute;
		%feature("autodoc", "	* /*! * \brief Computes mesh without geometry * \param aMesh - the mesh * \param aHelper - helper that must be used for adding elements to \aaMesh * etval bool - is a success * * The method is called if ( !aMesh->HasShapeToMesh() ) */

	:param aMesh:
	:type aMesh: SMESH_Mesh &
	:param aHelper:
	:type aHelper: SMESH_MesherHelper *
	:rtype: bool
") Compute;
		bool Compute (SMESH_Mesh & aMesh,SMESH_MesherHelper * aHelper);
		%feature("compactdefaultargs") Evaluate;
		%feature("autodoc", "	* /*! * \brief evaluates size of prospective mesh on a shape * \param aMesh - the mesh * \param aShape - the shape * \param aNbElems - prospective number of elements by types * etval bool - is a success */

	:param aMesh:
	:type aMesh: SMESH_Mesh &
	:param aShape:
	:type aShape: TopoDS_Shape &
	:param aResMap:
	:type aResMap: MapShapeNbElems &
	:rtype: bool
") Evaluate;
		bool Evaluate (SMESH_Mesh & aMesh,const TopoDS_Shape & aShape,MapShapeNbElems & aResMap);
		%feature("compactdefaultargs") GetUsedHypothesis;
		%feature("autodoc", "	* /*! * \brief Returns a list of compatible hypotheses used to mesh a shape * \param aMesh - the mesh * \param aShape - the shape * \param ignoreAuxiliary - do not include auxiliary hypotheses in the list * etval const std::list <const SMESHDS_Hypothesis*> - hypotheses list * * List the hypothesis used by the algorithm associated to the shape. * Hypothesis associated to father shape -are- taken into account (see * GetAppliedHypothesis). Relevant hypothesis have a name (type) listed in * the algorithm. This method could be surcharged by specific algorithms, in * case of several hypothesis simultaneously applicable. */

	:param aMesh:
	:type aMesh: SMESH_Mesh &
	:param aShape:
	:type aShape: TopoDS_Shape &
	:param ignoreAuxiliary: default value is true
	:type ignoreAuxiliary: bool
	:rtype: std::list< SMESHDS_Hypothesis *>
") GetUsedHypothesis;
		const std::list<const SMESHDS_Hypothesis *> & GetUsedHypothesis (SMESH_Mesh & aMesh,const TopoDS_Shape & aShape,const bool ignoreAuxiliary = true);
		%feature("compactdefaultargs") GetAppliedHypothesis;
		%feature("autodoc", "	* /*! * \brief Returns a list of compatible hypotheses assigned to a shape in a mesh * \param aMesh - the mesh * \param aShape - the shape * \param ignoreAuxiliary - do not include auxiliary hypotheses in the list * etval const std::list <const SMESHDS_Hypothesis*> - hypotheses list * * List the relevant hypothesis associated to the shape. Relevant hypothesis * have a name (type) listed in the algorithm. Hypothesis associated to * father shape -are not- taken into account (see GetUsedHypothesis) */

	:param aMesh:
	:type aMesh: SMESH_Mesh &
	:param aShape:
	:type aShape: TopoDS_Shape &
	:param ignoreAuxiliary: default value is true
	:type ignoreAuxiliary: bool
	:rtype: std::list< SMESHDS_Hypothesis *>
") GetAppliedHypothesis;
		const std::list<const SMESHDS_Hypothesis *> & GetAppliedHypothesis (SMESH_Mesh & aMesh,const TopoDS_Shape & aShape,const bool ignoreAuxiliary = true);
		%feature("compactdefaultargs") InitCompatibleHypoFilter;
		%feature("autodoc", "	* /*! * \brief Make the filter recognize only compatible hypotheses * \param theFilter - the filter to initialize * \param ignoreAuxiliary - make filter ignore compatible auxiliary hypotheses * etval bool - true if the algo has compatible hypotheses */

	:param theFilter:
	:type theFilter: SMESH_HypoFilter &
	:param ignoreAuxiliary:
	:type ignoreAuxiliary: bool
	:rtype: bool
") InitCompatibleHypoFilter;
		bool InitCompatibleHypoFilter (SMESH_HypoFilter & theFilter,const bool ignoreAuxiliary);
		%feature("compactdefaultargs") SetParametersByMesh;
		%feature("autodoc", "	* /*! * \brief Just return false as the algorithm does not hold parameters values */

	:param theMesh:
	:type theMesh: SMESH_Mesh *
	:param theShape:
	:type theShape: TopoDS_Shape &
	:rtype: bool
") SetParametersByMesh;
		bool SetParametersByMesh (const SMESH_Mesh * theMesh,const TopoDS_Shape & theShape);
		%feature("compactdefaultargs") SetParametersByDefaults;
		%feature("autodoc", "	:param dflts:
	:type dflts: SMESH_0D_Algo::TDefaults &
	:param theMesh: default value is 0
	:type theMesh: SMESH_Mesh *
	:rtype: bool
") SetParametersByDefaults;
		bool SetParametersByDefaults (const SMESH_0D_Algo::TDefaults & dflts,const SMESH_Mesh * theMesh = 0);
		%feature("compactdefaultargs") GetComputeError;
		%feature("autodoc", "	* /*! * \brief return compute error */

	:rtype: SMESH_ComputeErrorPtr
") GetComputeError;
		SMESH_ComputeErrorPtr GetComputeError ();
		%feature("compactdefaultargs") InitComputeError;
		%feature("autodoc", "	* /*! * \brief initialize compute error */

	:rtype: None
") InitComputeError;
		void InitComputeError ();
		%feature("compactdefaultargs") OnlyUnaryInput;
		%feature("autodoc", "	:rtype: bool
") OnlyUnaryInput;
		bool OnlyUnaryInput ();
		%feature("compactdefaultargs") NeedDescretBoundary;
		%feature("autodoc", "	:rtype: bool
") NeedDescretBoundary;
		bool NeedDescretBoundary ();
		%feature("compactdefaultargs") NeedShape;
		%feature("autodoc", "	:rtype: bool
") NeedShape;
		bool NeedShape ();
		%feature("compactdefaultargs") SupportSubmeshes;
		%feature("autodoc", "	:rtype: bool
") SupportSubmeshes;
		bool SupportSubmeshes ();
		%feature("compactdefaultargs") SetEventListener;
		%feature("autodoc", "	* /*! * \brief Sets event listener to submeshes if necessary * \param subMesh - submesh where algo is set * * This method is called when a submesh gets HYP_OK algo_state. * After being set, event listener is notified on each event of a submesh. * By default non listener is set */

	:param subMesh:
	:type subMesh: SMESH_subMesh *
	:rtype: None
") SetEventListener;
		void SetEventListener (SMESH_subMesh * subMesh);
		%feature("compactdefaultargs") SubmeshRestored;
		%feature("autodoc", "	* /*! * \brief Allow algo to do something after persistent restoration * \param subMesh - restored submesh * * This method is called only if a submesh has HYP_OK algo_state. */

	:param subMesh:
	:type subMesh: SMESH_subMesh *
	:rtype: None
") SubmeshRestored;
		void SubmeshRestored (SMESH_subMesh * subMesh);
		%feature("compactdefaultargs") GetNodeParamOnEdge;
		%feature("autodoc", "	* /*! * \brief Fill vector of node parameters on geometrical edge, including vertex nodes * \param theMesh - The mesh containing nodes * \param theEdge - The geometrical edge of interest * \param theParams - The resulting vector of sorted node parameters * etval bool - false if not all parameters are OK */

	:param theMesh:
	:type theMesh: SMESHDS_Mesh *
	:param theEdge:
	:type theEdge: TopoDS_Edge &
	:param theParams:
	:type theParams: std::vector< double> &
	:rtype: bool
") GetNodeParamOnEdge;
		static bool GetNodeParamOnEdge (const SMESHDS_Mesh * theMesh,const TopoDS_Edge & theEdge,std::vector< double> & theParams);
		%feature("compactdefaultargs") GetSortedNodesOnEdge;
		%feature("autodoc", "	* /*! * \brief Fill map of node parameter on geometrical edge to node it-self * \param theMesh - The mesh containing nodes * \param theEdge - The geometrical edge of interest * \param theNodes - The resulting map * \param ignoreMediumNodes - to store medium nodes of quadratic elements or not * etval bool - false if not all parameters are OK */

	:param theMesh:
	:type theMesh: SMESHDS_Mesh *
	:param theEdge:
	:type theEdge: TopoDS_Edge &
	:param ignoreMediumNodes:
	:type ignoreMediumNodes: bool
	:param theNodes:
	:type theNodes: std::map< double,  SMDS_MeshNode *> &
	:rtype: bool
") GetSortedNodesOnEdge;
		static bool GetSortedNodesOnEdge (const SMESHDS_Mesh * theMesh,const TopoDS_Edge & theEdge,const bool ignoreMediumNodes,std::map< double, const SMDS_MeshNode *> & theNodes);
		%feature("compactdefaultargs") IsReversedSubMesh;
		%feature("autodoc", "	* /*! * \brief Find out elements orientation on a geometrical face * \param theFace - The face correctly oriented in the shape being meshed * \param theMeshDS - The mesh data structure * etval bool - true if the face normal and the normal of first element * in the correspoding submesh point in different directions */

	:param theFace:
	:type theFace: TopoDS_Face &
	:param theMeshDS:
	:type theMeshDS: SMESHDS_Mesh *
	:rtype: bool
") IsReversedSubMesh;
		static bool IsReversedSubMesh (const TopoDS_Face & theFace,SMESHDS_Mesh * theMeshDS);
		%feature("compactdefaultargs") EdgeLength;
		%feature("autodoc", "	* /*! * \brief Compute length of an edge * \param E - the edge * etval double - the length */

	:param E:
	:type E: TopoDS_Edge &
	:rtype: double
") EdgeLength;
		static double EdgeLength (const TopoDS_Edge & E);
		%feature("compactdefaultargs") Continuity;
		%feature("autodoc", "	* /*! * \brief Return continuity of two edges * \param E1 - the 1st edge * \param E2 - the 2nd edge * etval GeomAbs_Shape - regularity at the junction between E1 and E2 */

	:param E1:
	:type E1: TopoDS_Edge &
	:param E2:
	:type E2: TopoDS_Edge &
	:rtype: GeomAbs_Shape
") Continuity;
		static GeomAbs_Shape Continuity (const TopoDS_Edge & E1,const TopoDS_Edge & E2);
		%feature("compactdefaultargs") IsContinuous;
		%feature("autodoc", "	* /*! * \brief Return true if an edge can be considered as a continuation of another */

	:param E1:
	:type E1: TopoDS_Edge &
	:param E2:
	:type E2: TopoDS_Edge &
	:rtype: bool
") IsContinuous;
		static bool IsContinuous (const TopoDS_Edge & E1,const TopoDS_Edge & E2);
		%feature("compactdefaultargs") VertexNode;
		%feature("autodoc", "	* /*! * \brief Return the node built on a vertex * \param V - the vertex * \param meshDS - mesh * etval const SMDS_MeshNode* - found node or NULL */

	:param V:
	:type V: TopoDS_Vertex &
	:param meshDS:
	:type meshDS: SMESHDS_Mesh *
	:rtype: SMDS_MeshNode *
") VertexNode;
		static const SMDS_MeshNode * VertexNode (const TopoDS_Vertex & V,const SMESHDS_Mesh * meshDS);
};


%nodefaultctor SMESH_HypoFilter;
class SMESH_HypoFilter : public SMESH_HypoPredicate {
/* public enums */
enum Logical {
	AND = 0,
	AND_NOT = 1,
	OR = 2,
	OR_NOT = 3,
};

enum Comparison {
	EQUAL = 0,
	NOT_EQUAL = 1,
	MORE = 2,
	LESS = 3,
};

/* end public enums declaration */

	public:
		%feature("compactdefaultargs") SMESH_HypoFilter;
		%feature("autodoc", "	:rtype: None
") SMESH_HypoFilter;
		 SMESH_HypoFilter ();
		%feature("compactdefaultargs") SMESH_HypoFilter;
		%feature("autodoc", "	:param aPredicate:
	:type aPredicate: SMESH_HypoPredicate *
	:param notNagate: default value is true
	:type notNagate: bool
	:rtype: None
") SMESH_HypoFilter;
		 SMESH_HypoFilter (SMESH_HypoPredicate * aPredicate,bool notNagate = true);
		%feature("compactdefaultargs") Init;
		%feature("autodoc", "	:param aPredicate:
	:type aPredicate: SMESH_HypoPredicate *
	:param notNagate: default value is true
	:type notNagate: bool
	:rtype: SMESH_HypoFilter
") Init;
		SMESH_HypoFilter & Init (SMESH_HypoPredicate * aPredicate,bool notNagate = true);
		%feature("compactdefaultargs") And;
		%feature("autodoc", "	:param aPredicate:
	:type aPredicate: SMESH_HypoPredicate *
	:rtype: SMESH_HypoFilter
") And;
		SMESH_HypoFilter & And (SMESH_HypoPredicate * aPredicate);
		%feature("compactdefaultargs") AndNot;
		%feature("autodoc", "	:param aPredicate:
	:type aPredicate: SMESH_HypoPredicate *
	:rtype: SMESH_HypoFilter
") AndNot;
		SMESH_HypoFilter & AndNot (SMESH_HypoPredicate * aPredicate);
		%feature("compactdefaultargs") Or;
		%feature("autodoc", "	:param aPredicate:
	:type aPredicate: SMESH_HypoPredicate *
	:rtype: SMESH_HypoFilter
") Or;
		SMESH_HypoFilter & Or (SMESH_HypoPredicate * aPredicate);
		%feature("compactdefaultargs") OrNot;
		%feature("autodoc", "	:param aPredicate:
	:type aPredicate: SMESH_HypoPredicate *
	:rtype: SMESH_HypoFilter
") OrNot;
		SMESH_HypoFilter & OrNot (SMESH_HypoPredicate * aPredicate);
		%feature("compactdefaultargs") IsAlgo;
		%feature("autodoc", "	:rtype: SMESH_HypoPredicate *
") IsAlgo;
		static SMESH_HypoPredicate * IsAlgo ();
		%feature("compactdefaultargs") IsAuxiliary;
		%feature("autodoc", "	:rtype: SMESH_HypoPredicate *
") IsAuxiliary;
		static SMESH_HypoPredicate * IsAuxiliary ();
		%feature("compactdefaultargs") IsApplicableTo;
		%feature("autodoc", "	:param theShape:
	:type theShape: TopoDS_Shape &
	:rtype: SMESH_HypoPredicate *
") IsApplicableTo;
		static SMESH_HypoPredicate * IsApplicableTo (const TopoDS_Shape & theShape);
		%feature("compactdefaultargs") IsAssignedTo;
		%feature("autodoc", "	:param theShape:
	:type theShape: TopoDS_Shape &
	:rtype: SMESH_HypoPredicate *
") IsAssignedTo;
		static SMESH_HypoPredicate * IsAssignedTo (const TopoDS_Shape & theShape);
		%feature("compactdefaultargs") Is;
		%feature("autodoc", "	:param theHypo:
	:type theHypo: SMESH_Hypothesis *
	:rtype: SMESH_HypoPredicate *
") Is;
		static SMESH_HypoPredicate * Is (const SMESH_Hypothesis * theHypo);
		%feature("compactdefaultargs") IsGlobal;
		%feature("autodoc", "	:param theMainShape:
	:type theMainShape: TopoDS_Shape &
	:rtype: SMESH_HypoPredicate *
") IsGlobal;
		static SMESH_HypoPredicate * IsGlobal (const TopoDS_Shape & theMainShape);
		%feature("compactdefaultargs") IsMoreLocalThan;
		%feature("autodoc", "	:param theShape:
	:type theShape: TopoDS_Shape &
	:rtype: SMESH_HypoPredicate *
") IsMoreLocalThan;
		static SMESH_HypoPredicate * IsMoreLocalThan (const TopoDS_Shape & theShape);
		%feature("compactdefaultargs") HasName;
		%feature("autodoc", "	:param theName:
	:type theName: std::string &
	:rtype: SMESH_HypoPredicate *
") HasName;
		static SMESH_HypoPredicate * HasName (const std::string & theName);
		%feature("compactdefaultargs") HasDim;
		%feature("autodoc", "	:param theDim:
	:type theDim: int
	:rtype: SMESH_HypoPredicate *
") HasDim;
		static SMESH_HypoPredicate * HasDim (const int theDim);
		%feature("compactdefaultargs") HasType;
		%feature("autodoc", "	:param theHypType:
	:type theHypType: int
	:rtype: SMESH_HypoPredicate *
") HasType;
		static SMESH_HypoPredicate * HasType (const int theHypType);
		%feature("compactdefaultargs") IsOk;
		%feature("autodoc", "	* /*! * \brief check aHyp or/and aShape it is assigned to */

	:param aHyp:
	:type aHyp: SMESH_Hypothesis *
	:param aShape:
	:type aShape: TopoDS_Shape &
	:rtype: bool
") IsOk;
		bool IsOk (const SMESH_Hypothesis * aHyp,const TopoDS_Shape & aShape);
		%feature("compactdefaultargs") IsAny;
		%feature("autodoc", "	* /*! * \brief return true if contains no predicates */

	:rtype: bool
") IsAny;
		bool IsAny ();
};


%nodefaultctor SMESH_OctreeNode;
class SMESH_OctreeNode : public SMESH_Octree {
	public:
		%feature("compactdefaultargs") SMESH_OctreeNode;
		%feature("autodoc", "	:param theNodes:
	:type theNodes: std::set< SMDS_MeshNode *> &
	:param maxLevel: default value is -1
	:type maxLevel: int
	:param maxNbNodes: default value is 5
	:type maxNbNodes: int
	:param minBoxSize: default value is 0
	:type minBoxSize: double
	:rtype: None
") SMESH_OctreeNode;
		 SMESH_OctreeNode (const std::set<const SMDS_MeshNode *> & theNodes,const int maxLevel = -1,const int maxNbNodes = 5,const double minBoxSize = 0);
		%feature("compactdefaultargs") isInside;
		%feature("autodoc", "	:param Node:
	:type Node: SMDS_MeshNode *
	:param precision: default value is 0
	:type precision: double
	:rtype: bool
") isInside;
		const bool isInside (const SMDS_MeshNode * Node,const double precision = 0);
		%feature("compactdefaultargs") NodesAround;
		%feature("autodoc", "	:param Node:
	:type Node: SMDS_MeshNode *
	:param Result:
	:type Result: std::list< SMDS_MeshNode *> *
	:param precision: default value is 0
	:type precision: double
	:rtype: None
") NodesAround;
		void NodesAround (const SMDS_MeshNode * Node,std::list<const SMDS_MeshNode *> * Result,const double precision = 0);
		%feature("compactdefaultargs") NodesAround;
		%feature("autodoc", "	:param Node:
	:type Node: SMDS_MeshNode *
	:param dist2Nodes:
	:type dist2Nodes: std::map<double,  SMDS_MeshNode *> &
	:param precision:
	:type precision: double
	:rtype: bool
") NodesAround;
		bool NodesAround (const SMDS_MeshNode * Node,std::map<double, const SMDS_MeshNode *> & dist2Nodes,double precision);
		%feature("compactdefaultargs") FindCoincidentNodes;
		%feature("autodoc", "	:param nodes:
	:type nodes: std::set< SMDS_MeshNode *> *
	:param theTolerance:
	:type theTolerance: double
	:param theGroupsOfNodes:
	:type theGroupsOfNodes: std::list< std::list<  SMDS_MeshNode *> > *
	:rtype: None
") FindCoincidentNodes;
		void FindCoincidentNodes (std::set<const SMDS_MeshNode *> * nodes,const double theTolerance,std::list< std::list< const SMDS_MeshNode *> > * theGroupsOfNodes);
		%feature("compactdefaultargs") FindCoincidentNodes;
		%feature("autodoc", "	:param nodes:
	:type nodes: std::set< SMDS_MeshNode *> &
	:param theGroupsOfNodes:
	:type theGroupsOfNodes: std::list< std::list<  SMDS_MeshNode *> > *
	:param theTolerance: default value is 0.00001
	:type theTolerance: double
	:param maxLevel: default value is -1
	:type maxLevel: int
	:param maxNbNodes: default value is 5
	:type maxNbNodes: int
	:rtype: void
") FindCoincidentNodes;
		static void FindCoincidentNodes (std::set<const SMDS_MeshNode *> & nodes,std::list< std::list< const SMDS_MeshNode *> > * theGroupsOfNodes,const double theTolerance = 0.00001,const int maxLevel = -1,const int maxNbNodes = 5);
		%feature("compactdefaultargs") UpdateByMoveNode;
		%feature("autodoc", "	* /*! * \brief Update data according to node movement */

	:param node:
	:type node: SMDS_MeshNode *
	:param toPnt:
	:type toPnt: gp_Pnt
	:rtype: None
") UpdateByMoveNode;
		void UpdateByMoveNode (const SMDS_MeshNode * node,const gp_Pnt & toPnt);
		%feature("compactdefaultargs") GetChildrenIterator;
		%feature("autodoc", "	* /*! * \brief Return iterator over children */

	:rtype: SMESH_OctreeNodeIteratorPtr
") GetChildrenIterator;
		SMESH_OctreeNodeIteratorPtr GetChildrenIterator ();
		%feature("compactdefaultargs") GetNodeIterator;
		%feature("autodoc", "	* /*! * \brief Return nodes iterator */

	:rtype: SMDS_NodeIteratorPtr
") GetNodeIterator;
		SMDS_NodeIteratorPtr GetNodeIterator ();
		%feature("compactdefaultargs") NbNodes;
		%feature("autodoc", "	* /*! * \brief Return nb nodes in a tree */

	:rtype: int
") NbNodes;
		int NbNodes ();
};


%nodefaultctor SMESH_0D_Algo;
class SMESH_0D_Algo : public SMESH_Algo {
	public:
		%feature("compactdefaultargs") SMESH_0D_Algo;
		%feature("autodoc", "	:param hypId:
	:type hypId: int
	:param studyId:
	:type studyId: int
	:param gen:
	:type gen: SMESH_Gen *
	:rtype: None
") SMESH_0D_Algo;
		 SMESH_0D_Algo (int hypId,int studyId,SMESH_Gen * gen);
};


%nodefaultctor SMESH_1D_Algo;
class SMESH_1D_Algo : public SMESH_Algo {
	public:
		%feature("compactdefaultargs") SMESH_1D_Algo;
		%feature("autodoc", "	:param hypId:
	:type hypId: int
	:param studyId:
	:type studyId: int
	:param gen:
	:type gen: SMESH_Gen *
	:rtype: None
") SMESH_1D_Algo;
		 SMESH_1D_Algo (int hypId,int studyId,SMESH_Gen * gen);
};


%nodefaultctor SMESH_2D_Algo;
class SMESH_2D_Algo : public SMESH_Algo {
	public:
		%feature("compactdefaultargs") SMESH_2D_Algo;
		%feature("autodoc", "	:param hypId:
	:type hypId: int
	:param studyId:
	:type studyId: int
	:param gen:
	:type gen: SMESH_Gen *
	:rtype: None
") SMESH_2D_Algo;
		 SMESH_2D_Algo (int hypId,int studyId,SMESH_Gen * gen);
		%feature("compactdefaultargs") NumberOfWires;
		%feature("autodoc", "	:param S:
	:type S: TopoDS_Shape &
	:rtype: int
") NumberOfWires;
		int NumberOfWires (const TopoDS_Shape & S);
		%feature("compactdefaultargs") NumberOfPoints;
		%feature("autodoc", "	:param aMesh:
	:type aMesh: SMESH_Mesh &
	:param W:
	:type W: TopoDS_Wire &
	:rtype: int
") NumberOfPoints;
		int NumberOfPoints (SMESH_Mesh & aMesh,const TopoDS_Wire & W);
};


%nodefaultctor SMESH_3D_Algo;
class SMESH_3D_Algo : public SMESH_Algo {
	public:
		%feature("compactdefaultargs") SMESH_3D_Algo;
		%feature("autodoc", "	:param hypId:
	:type hypId: int
	:param studyId:
	:type studyId: int
	:param gen:
	:type gen: SMESH_Gen *
	:rtype: None
") SMESH_3D_Algo;
		 SMESH_3D_Algo (int hypId,int studyId,SMESH_Gen * gen);
};


