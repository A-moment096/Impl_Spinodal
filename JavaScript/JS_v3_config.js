import fs from "node:fs";

class ConfigValueReader {
  static number(value, keyPath) {
    if (typeof value !== "number" || Number.isNaN(value)) {
      throw new Error(`Config '${keyPath}' must be a valid number.`);
    }
    return value;
  }

  static string(value, keyPath) {
    if (typeof value !== "string" || value.trim() === "") {
      throw new Error(`Config '${keyPath}' must be a non-empty string.`);
    }
    return value;
  }
}

class MeshConfig {
  constructor(Nx, Ny, dx, dy, con_init, dcon, boundaryConditions) {
    this.Nx = Nx;
    this.Ny = Ny;
    this.dx = dx;
    this.dy = dy;
    this.con_init = con_init;
    this.dcon = dcon;
    this.boundaryConditions = boundaryConditions;
  }

  static fromParsed(parsed) {
    if (!Array.isArray(parsed?.meshes) || parsed.meshes.length === 0) {
      throw new Error("Config 'meshes' must be a non-empty array.");
    }

    const targetMesh =
      parsed.meshes.find((mesh) => mesh?.type === "Cahn-Hilliard") ??
      parsed.meshes[0];

    const size = targetMesh?.size;
    if (!Array.isArray(size) || size.length !== 2) {
      throw new Error("Config 'meshes[...].size' must be [Nx, Ny].");
    }

    const spacing = targetMesh?.spacing;
    if (!Array.isArray(spacing) || spacing.length !== 2) {
      throw new Error("Config 'meshes[...].spacing' must be [dx, dy].");
    }

    const boundaryConditions = targetMesh?.boundary_conditions;
    if (!Array.isArray(boundaryConditions) || boundaryConditions.length === 0) {
      throw new Error(
        "Config 'meshes[...].boundary_conditions' must be a non-empty array.",
      );
    }

    const initType = ConfigValueReader.string(
      targetMesh?.initial_geometry?.type,
      "meshes[...].initial_geometry.type",
    );
    if (initType !== "random") {
      throw new Error(
        `Unsupported initial geometry type '${initType}'. Only 'random' is supported.`,
      );
    }

    const minVal = ConfigValueReader.number(
      targetMesh?.initial_geometry?.parameters?.min_val,
      "meshes[...].initial_geometry.parameters.min_val",
    );
    const maxVal = ConfigValueReader.number(
      targetMesh?.initial_geometry?.parameters?.max_val,
      "meshes[...].initial_geometry.parameters.max_val",
    );
    if (maxVal < minVal) {
      throw new Error(
        "Config 'meshes[...].initial_geometry.parameters.max_val' must be >= min_val.",
      );
    }

    return new MeshConfig(
      ConfigValueReader.number(size[0], "meshes[...].size[0]"),
      ConfigValueReader.number(size[1], "meshes[...].size[1]"),
      ConfigValueReader.number(spacing[0], "meshes[...].spacing[0]"),
      ConfigValueReader.number(spacing[1], "meshes[...].spacing[1]"),
      0.5 * (minVal + maxVal),
      0.5 * (maxVal - minVal),
      boundaryConditions,
    );
  }
}

class EnergyConfig {
  constructor(A, kappa) {
    this.A = A;
    this.kappa = kappa;
  }

  static fromParsed(parsed) {
    const bulkType = ConfigValueReader.string(
      parsed?.energy?.bulk_energy?.type,
      "energy.bulk_energy.type",
    );
    if (bulkType !== "double_well") {
      throw new Error(
        `Unsupported bulk energy type '${bulkType}'. Only 'double_well' is supported.`,
      );
    }

    const interfacialType = ConfigValueReader.string(
      parsed?.energy?.interfacial_energy?.type,
      "energy.interfacial_energy.type",
    );
    if (interfacialType !== "gradient") {
      throw new Error(
        `Unsupported interfacial energy type '${interfacialType}'. Only 'gradient' is supported.`,
      );
    }

    return new EnergyConfig(
      ConfigValueReader.number(
        parsed?.energy?.bulk_energy?.parameters?.A,
        "energy.bulk_energy.parameters.A",
      ),
      ConfigValueReader.number(
        parsed?.energy?.interfacial_energy?.parameters?.kappa,
        "energy.interfacial_energy.parameters.kappa",
      ),
    );
  }
}

class SolverConfig {
  constructor(Nstep, dt, M) {
    this.Nstep = Nstep;
    this.dt = dt;
    this.M = M;
  }

  static fromParsed(parsed) {
    const solverList = parsed?.solvers?.solver;
    if (!Array.isArray(solverList) || solverList.length === 0) {
      throw new Error("Config 'solvers.solver' must be a non-empty array.");
    }
    const chSolver = solverList.find(
      (solver) => solver?.type === "Cahn-Hilliard",
    );
    if (!chSolver) {
      throw new Error(
        "Config 'solvers.solver' must include a Cahn-Hilliard solver.",
      );
    }

    return new SolverConfig(
      ConfigValueReader.number(parsed?.solvers?.Nstep, "solvers.Nstep"),
      ConfigValueReader.number(parsed?.solvers?.dt, "solvers.dt"),
      ConfigValueReader.number(chSolver?.M, "solvers.solver[...].M"),
    );
  }
}

class OutputConfig {
  constructor(folder_name, writeEvery, logEvery, nanCheckEvery) {
    this.folder_name = folder_name;
    this.writeEvery = writeEvery;
    this.logEvery = logEvery;
    this.nanCheckEvery = nanCheckEvery;
  }

  static fromParsed(parsed) {
    const nanCheckEvery = parsed?.output?.nanCheckEvery;
    if (
      nanCheckEvery !== false &&
      (typeof nanCheckEvery !== "number" ||
        Number.isNaN(nanCheckEvery) ||
        !Number.isInteger(nanCheckEvery) ||
        nanCheckEvery <= 0)
    ) {
      throw new Error(
        "Config 'output.nanCheckEvery' must be a positive integer or false.",
      );
    }

    return new OutputConfig(
      ConfigValueReader.string(parsed?.output?.folder, "output.folder"),
      ConfigValueReader.number(parsed?.output?.writeEvery, "output.writeEvery"),
      ConfigValueReader.number(parsed?.output?.logEvery, "output.logEvery"),
      nanCheckEvery,
    );
  }
}

class SimulationConfig {
  constructor(mesh, energy, solver, output) {
    this.mesh = mesh;
    this.energy = energy;
    this.solver = solver;
    this.output = output;
  }

  static fromFile(configPath) {
    let parsed;
    try {
      const fileText = fs.readFileSync(configPath, "utf-8");
      parsed = JSON.parse(fileText);
    } catch (error) {
      throw new Error(
        `Failed to load config file '${configPath}': ${error.message}`,
      );
    }

    return new SimulationConfig(
      MeshConfig.fromParsed(parsed),
      EnergyConfig.fromParsed(parsed),
      SolverConfig.fromParsed(parsed),
      OutputConfig.fromParsed(parsed),
    );
  }
}

export { SimulationConfig };
