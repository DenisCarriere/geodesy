import path from 'path'
import resolve from 'rollup-plugin-node-resolve'
import commonjs from 'rollup-plugin-commonjs'
import builtins from 'rollup-plugin-node-builtins'
import globals from 'rollup-plugin-node-globals'
import json from 'rollup-plugin-json'
import cleanup from 'rollup-plugin-cleanup'
const pckg = require('./package.json')

export default {
  entry: pckg['jsnext:main'],
  dest: pckg['main'],
  format: 'cjs',
  plugins: [
    cleanup(),
    json(),
    resolve(),
    commonjs(),
    globals(),
    builtins()
  ],
  useStrict: false,
}