/** Which fit produced a dab. Drives both its colour and its visibility toggle. */
export const DabKind = {
  CLOTHOID: 0,
  BEZIER  : 1,
  INPUT   : 2,
};

const STRIDE = 4;

/**
 * A flat x/y/radius/kind buffer of brush dabs.
 *
 * Flat rather than an array of objects because a single fast stroke lays down tens of
 * thousands of dabs and the demo never removes one individually.
 */
export class DabList {
  data: number[] = [];

  get count() {
    return this.data.length / STRIDE;
  }

  push(x: number, y: number, radius: number, kind: number) {
    this.data.push(x, y, radius, kind);
  }

  clear() {
    this.data.length = 0;
  }

  forEach(cb: (x: number, y: number, radius: number, kind: number) => void) {
    const d = this.data;

    for (let i = 0; i < d.length; i += STRIDE) {
      cb(d[i], d[i + 1], d[i + 2], d[i + 3]);
    }
  }
}
