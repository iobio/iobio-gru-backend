import assert from 'node:assert/strict';
import test from 'node:test';
import { applyConfigEnvOverrides, parseEnvValue } from '../src/config.js';

test('maps public IOBIO environment namespaces to snake_case config keys', () => {
  const config = {
    gene: { path_prefix: '/gene', show_intro: false },
    backend: { path_prefix: '/' },
  };

  applyConfigEnvOverrides(config, {
    IOBIO_GENE_DEFAULT_MODE: 'simple',
    IOBIO_GENE_SHOW_INTRO: 'true',
    IOBIO_GENE_SHOW_FILES_BUTTON: 'false',
    IOBIO_GENE_SITE_NAME: 'Nebula Genomics',
    IOBIO_GENE_CUSTOM_SETTING: 'published',
    IOBIO_BACKEND_PATH_PREFIX: '/gru/api/v1',
    IOBIO_SECRET: 'not public config',
    OTHER_GENE_SHOW_INTRO: 'false',
  });

  assert.deepEqual(config, {
    gene: {
      path_prefix: '/gene',
      show_intro: true,
      default_mode: 'simple',
      show_files_button: false,
      site_name: 'Nebula Genomics',
      custom_setting: 'published',
    },
    backend: { path_prefix: '/gru/api/v1' },
  });
});

test('creates a supported config section when it is absent', () => {
  const config = {};
  applyConfigEnvOverrides(config, { IOBIO_BAM_PATH_PREFIX: '/bam' });
  assert.deepEqual(config, { bam: { path_prefix: '/bam' } });
});

test('only exact lowercase boolean strings are converted', () => {
  assert.equal(parseEnvValue('true'), true);
  assert.equal(parseEnvValue('false'), false);
  assert.equal(parseEnvValue('TRUE'), 'TRUE');
  assert.equal(parseEnvValue('123'), '123');
  assert.equal(parseEnvValue(''), '');
});
