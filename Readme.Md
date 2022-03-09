Polynomial clothoid-based spline solver that's
optmized for brush strokes.

WIP.

General ideas:
* Try out shortest path decomposition paper?
 - Would be fairly easy to do with a bspline curvature function
 - No need to have seperate arc/clothoid/line primitives in that
   case, you'd just use a degree 1 (linear) bspline.

# Setup

```
git clone https://github.com/joeedh/research-framework.git
cd research-framework
npm run-script setup

node serv.js
```


