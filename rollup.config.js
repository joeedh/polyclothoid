function isExternal(id) {
  return id.search(/path\.ux/) >= 0 || id.search(/pathux/) >= 0;
}

export default {
  input: 'scripts/stroker/stroke.js',
  output: {
    file: 'build/stroker.js',
    format: 'es',
  },
  external : id => (console.log(id, isExternal(id)), isExternal(id)),
  inlineDynamicImports : true
};
