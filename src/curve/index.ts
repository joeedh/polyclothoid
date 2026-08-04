export { Curve, type Canvas2DLike, type CanvasPaint, type CurveEdge, type CurveVertex } from "./curve.js";
export {
  type CurveSolver,
  type Pairing,
  type SolvableEdge,
  type SolvableMesh,
  type SolvableVertex,
} from "./mesh_types.js";
export { chains, walk, type Chain } from "./topology.js";

export {
  defaultLevel,
  findPairing,
  levelFor,
  maxLevel,
  pairingLevel,
  sharing,
  vertexPartitions,
  type Sharing,
} from "./pairing.js";

export {
  RATE_WINDOW,
  diagnosticThresholds,
  emptySolveReport,
  geometricRate,
  stabilityThresholds,
  type Diagnostic,
  type DiagnosticAction,
  type DiagnosticCondition,
  type DiagnosticSeverity,
  type SolveReport,
  type SolveTrace,
  type TraceStep,
} from "./diagnostics.js";

export {
  Clothoid,
  ClothoidSolver,
  KORDER,
  KSCALE,
  KTH,
  KOFFX,
  KOFFY,
  KARCSCALE,
  KTOT,
  defaultClothoidSolverOptions,
  type ClothoidSolverOptions,
} from "./clothoid.js";

export {
  activeProfile,
  setCurvatureProfile,
  piecewiseLinear,
  circleArc,
  sPowerProfile,
  bernsteinCurvature,
  profileByName,
  registerProfile,
  type CurvatureProfile,
} from "./profile.js";

export {
  MAX_SPOWER_ORDER,
  differentiateSPower,
  differentiationMatrix,
  endpointTaylor,
  evalSPower,
  hik,
  integralWeights,
  integrateSPower,
  massMatrix,
  pairsToTaylor,
  reverseCurvature,
  reverseSPower,
  sPowerCurvature2,
  sPowerDerivative,
  sPowerIntegral,
  sPowerLength,
  sPowerOrder,
  stiffnessMatrix,
  taylorToPairs,
} from "./spower.js";

export {
  MIN_CANONICAL_CHORD,
  defaultQuadratureFor,
  defaultQuadratureOptions,
  defaultSPowerQuadrature,
  integrateProfile,
  integrateProfileJacobian,
  totalTurning,
  transformJacobian,
  type QuadratureOptions,
  type TransformJacobian,
} from "./quadrature.js";

export {
  blockLength,
  congruence,
  continuationEntry,
  edgeCoefficients,
  edgeDOFLength,
  edgeEnergy,
  edgeHessian,
  edgePullback,
  edgeTransform,
  edgeTurningRow,
  matrixTimesVector,
  multiply,
  referenceLength,
  rowTimesMatrix,
  sPowerDOF,
  transformEntry,
  type EdgeFrame,
  type ProfileDOF,
} from "./blocks.js";

export {
  SAMPLE_COUNT,
  sampleDOF,
  sampleEnergy,
  sampleMatrix,
  sampleTransform,
  sampledDOF,
  trapezoidWeights,
} from "./samples.js";

export { SPOWER_ORDER, SPowerClothoid } from "./spower_clothoid.js";

export {
  SPowerSolver,
  defaultSPowerSolverOptions,
  type ChainRun,
  type Fault,
  type SPowerSolverOptions,
  type SPowerSolverReport,
} from "./spower_solver.js";

export { CubicBezier, BezierSolver, cubic } from "./bezier.js";
export { BSpline, BSplinePoint, BSplineSolver, BSplineRecalc, dbasis } from "./bspline.js";
