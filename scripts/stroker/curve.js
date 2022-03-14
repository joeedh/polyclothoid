export class Curve {
  evaluate(s) {
    throw new Error("implement me!");
  }

  derivative(s) {
    throw new Error("implement me!");
  }

  derivative2(s) {
    throw new Error("implement me!");
  }

  curvature(s) {
    let dv1 = this.derivative(s);
    let dv2 = this.derivative2(s);

    return (dv1[0]*dv2[1] - dv1[1]*dv2[0]) / Math.pow(dv1.dot(dv1), 3.0/2.0);
  }

  distanceTo(s) {
    throw new Error("implement me!");
  }

  get length() {
    throw new Error("implement me!");
  }
}