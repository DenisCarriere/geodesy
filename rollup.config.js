import cleanup from 'rollup-plugin-cleanup'

export default {
  entry: 'index.js',
  plugins: [
    cleanup()
  ],
  targets: [
    {
      format: 'es',
      dest: 'dist/index.es6.js'
    },
    {
      format: 'cjs',
      dest: 'dist/index.js'
    }
  ],
  useStrict: false,
};