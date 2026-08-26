const PUBLIC_CONFIG_SECTIONS = new Set(['bam', 'gene', 'backend']);

export function applyConfigEnvOverrides(config, env = process.env) {
  for (const [name, rawValue] of Object.entries(env)) {
    const match = /^IOBIO_([A-Z0-9]+)_(.+)$/.exec(name);
    if (!match) {
      continue;
    }

    const section = match[1].toLowerCase();
    if (!PUBLIC_CONFIG_SECTIONS.has(section)) {
      continue;
    }

    const key = match[2].toLowerCase();
    config[section] = config[section] || {};
    config[section][key] = parseEnvValue(rawValue);
  }

  return config;
}

export function parseEnvValue(value) {
  if (value === 'true') {
    return true;
  }
  if (value === 'false') {
    return false;
  }
  return value;
}
