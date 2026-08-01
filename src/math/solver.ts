import { listRemove } from "./util.js";

/** A block of solver-owned scalars. In practice a `Float64Array` view over a `ks` array. */
export type ParamVector = Float64Array | number[];

/**
 * A single constraint in the {@link Solver}.
 *
 * A constraint reports a scalar error for the current parameter values. The solver
 * finite-differences that error to get a gradient and steps the parameters down it.
 *
 * **Returning `0.0` is meaningful and load-bearing.** It short-circuits twice: here in
 * {@link evaluate}, before any gradient is finite-differenced, and again in
 * {@link Solver.solveStep}, before any parameter is touched. A constraint that returns
 * `0.0` unconditionally is therefore invoked once per iteration purely for its side
 * effects, at the cost of one call rather than `2 x paramCount` evaluations.
 *
 * The clothoid solver's curvature constraint relies on exactly this to act as a direct
 * projection interleaved Gauss-Seidel style with the gradient descent of the other
 * constraints. Do not "fix" it into returning a residual.
 */
export class Constraint<P = unknown> {
  glst: Float64Array[] = [];
  wlst: Float64Array[] = [];

  df = 0.0005;
  threshold = 0.0001;

  /** Optional analytic gradient. When set, finite differencing is skipped entirely. */
  funcDv: ((params: P, glst: Float64Array[]) => void) | null = null;

  constructor(
    public name: string,
    public func: (params: P) => number,
    public klst: ParamVector[],
    public params: P,
    public k = 1.0
  ) {
    for (const ks of klst) {
      this.glst.push(new Float64Array(ks.length));

      const ws = new Float64Array(ks.length);
      ws.fill(1.0);
      this.wlst.push(ws);
    }
  }

  /** Hook for constraints that need to normalize their parameters after a step. */
  postSolve() {
    /* no-op */
  }

  evaluate(noDvs = false) {
    const r1 = this.func(this.params);

    if (this.funcDv) {
      this.funcDv(this.params, this.glst);
      return r1;
    }

    if (Math.abs(r1) < this.threshold) {
      return 0.0;
    }

    if (noDvs) {
      return r1;
    }

    const df = this.df;

    for (let i = 0; i < this.klst.length; i++) {
      const gs = this.glst[i];
      const ks = this.klst[i];

      for (let j = 0; j < ks.length; j++) {
        const orig = ks[j];

        ks[j] += df;
        const r2 = this.func(this.params);
        ks[j] = orig;

        gs[j] = (r2 - r1) / df;
      }
    }

    return r1;
  }
}

/**
 * Gauss-Seidel style solver: constraints are visited in order and each one steps the
 * parameters immediately, so later constraints in an iteration see earlier updates.
 */
export class Solver {
  constraints: Constraint<never>[] = [];

  gk = 0.99;
  simple = false;
  randCons = false;

  /** Stop once total constraint error falls below this. */
  threshold = 0.01;

  add<P>(con: Constraint<P>) {
    this.constraints.push(con as unknown as Constraint<never>);
  }

  remove<P>(con: Constraint<P>) {
    listRemove(this.constraints, con as unknown as Constraint<never>);
  }

  private pick(i: number) {
    if (this.randCons) {
      return this.constraints[~~(Math.random() * this.constraints.length * 0.99999)];
    }

    return this.constraints[i];
  }

  solveStep(gk = this.gk) {
    let err = 0.0;

    for (let ci = 0; ci < this.constraints.length; ci++) {
      const con = this.pick(ci);
      let r1 = con.evaluate();

      // Side-effect-only constraints report 0.0; skip before touching any parameter.
      if (r1 === 0.0) {
        continue;
      }

      err += Math.abs(r1);

      let totgs = 0.0;
      for (let i = 0; i < con.klst.length; i++) {
        const gs = con.glst[i];

        for (let j = 0; j < con.klst[i].length; j++) {
          totgs += gs[j] * gs[j];
        }
      }

      if (totgs === 0.0) {
        continue;
      }

      r1 /= totgs;

      for (let i = 0; i < con.klst.length; i++) {
        const ks = con.klst[i];
        const gs = con.glst[i];
        const ws = con.wlst[i];

        for (let j = 0; j < ks.length; j++) {
          ks[j] += -r1 * gs[j] * con.k * gk * ws[j];
        }
      }

      con.postSolve();
    }

    return err;
  }

  /** Fixed-magnitude step along the gradient rather than a Newton-ish scaled one. */
  solveStepSimple(gk = this.gk) {
    let err = 0.0;

    for (let ci = 0; ci < this.constraints.length; ci++) {
      const con = this.pick(ci);
      const r1 = con.evaluate();

      if (r1 === 0.0) {
        continue;
      }

      err += Math.abs(r1);

      let totgs = 0.0;
      for (let i = 0; i < con.klst.length; i++) {
        const gs = con.glst[i];

        for (let j = 0; j < con.klst[i].length; j++) {
          totgs += gs[j] * gs[j];
        }
      }

      if (totgs === 0.0) {
        continue;
      }

      totgs = 0.0001 / Math.sqrt(totgs);

      for (let i = 0; i < con.klst.length; i++) {
        const ks = con.klst[i];
        const gs = con.glst[i];
        const ws = con.wlst[i];

        for (let j = 0; j < ks.length; j++) {
          ks[j] += -totgs * gs[j] * con.k * gk * ws[j];
        }
      }

      con.postSolve();
    }

    return err;
  }

  solve(steps: number, gk = this.gk, printError = false) {
    let err = 0.0;

    for (let i = 0; i < steps; i++) {
      err = this.simple ? this.solveStepSimple(gk) : this.solveStep(gk);

      if (printError) {
        console.warn("average error:", (err / this.constraints.length).toFixed(4));
      }

      if (err < this.threshold / this.constraints.length) {
        break;
      }
    }

    return err;
  }
}
