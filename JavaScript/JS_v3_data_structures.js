import fs from "node:fs";

class BoundaryCondition {
  constructor(pos, type) {
    this.type = type;
    this.pos = pos;
  }
}

class PeriodicBC extends BoundaryCondition {
  constructor(pos) {
    super(pos, "periodic");
  }
  apply(mesh) {
    const stride = mesh.Nx + 2;
    if (this.pos === "north") {
      for (let idx = 0; idx < mesh.Nx; idx++) {
        mesh.data[(mesh.Ny + 1) * stride + idx + 1] =
          mesh.data[1 * stride + idx + 1];
      }
    }
    if (this.pos === "south") {
      for (let idx = 0; idx < mesh.Nx; idx++) {
        mesh.data[idx + 1] = mesh.data[mesh.Ny * stride + idx + 1];
      }
    }
    if (this.pos === "east") {
      for (let idx = 0; idx < mesh.Ny; idx++) {
        mesh.data[(idx + 1) * stride + (mesh.Nx + 1)] =
          mesh.data[(idx + 1) * stride + 1];
      }
    }
    if (this.pos === "west") {
      for (let idx = 0; idx < mesh.Ny; idx++) {
        mesh.data[(idx + 1) * stride] = mesh.data[(idx + 1) * stride + mesh.Nx];
      }
    }
  }
}

class FixedBC extends BoundaryCondition {
  constructor(pos, value) {
    super(pos, "fixed");
    this.value = value;
  }
  apply(mesh) {
    const stride = mesh.Nx + 2;
    if (this.pos === "north") {
      for (let idx = 0; idx < mesh.Nx; idx++) {
        mesh.data[(mesh.Ny + 1) * stride + idx + 1] = this.value;
      }
    }
    if (this.pos === "south") {
      for (let idx = 0; idx < mesh.Nx; idx++) {
        mesh.data[idx + 1] = this.value;
      }
    }
    if (this.pos === "east") {
      for (let idx = 0; idx < mesh.Ny; idx++) {
        mesh.data[(idx + 1) * stride + (mesh.Nx + 1)] = this.value;
      }
    }
    if (this.pos === "west") {
      for (let idx = 0; idx < mesh.Ny; idx++) {
        mesh.data[(idx + 1) * stride] = this.value;
      }
    }
  }
}

class BoundaryConditions {
  constructor(bcConfigs = undefined) {
    this.bcs = [];
    if (bcConfigs !== undefined) {
      this.loadFromConfigList(bcConfigs);
    }
  }

  createFromConfig(bcConfig) {
    const position = bcConfig?.position;
    const type = bcConfig?.type;

    if (!["north", "south", "east", "west"].includes(position)) {
      throw new Error(
        `Invalid boundary position '${position}'. Use north/south/east/west.`,
      );
    }

    if (type === "periodic") {
      return new PeriodicBC(position);
    }

    if (type === "fixed") {
      if (typeof bcConfig?.value !== "number" || Number.isNaN(bcConfig.value)) {
        throw new Error(
          `Boundary '${position}' with type 'fixed' requires numeric 'value'.`,
        );
      }
      return new FixedBC(position, bcConfig.value);
    }

    throw new Error(`Unsupported boundary type '${type}' at '${position}'.`);
  }

  loadFromConfigList(bcConfigs) {
    if (!Array.isArray(bcConfigs)) {
      throw new Error("Boundary condition config must be an array.");
    }
    bcConfigs.forEach((bcConfig) => this.addBC(this.createFromConfig(bcConfig)));
    return this;
  }

  addBC(bc) {
    this.bcs.push(bc);
  }
  applyAll(mesh) {
    this.bcs.forEach((bc) => bc.apply(mesh));
  }
}

class InitGeometry {
  constructor(type, params) {
    this.type = type;
    this.params = params;
  }
}

class RandomInit extends InitGeometry {
  constructor(baseValue, minVar, maxVar) {
    super("random", { baseValue, minVar, maxVar });
  }
  apply(mesh) {
    const stride = mesh.Nx + 2;
    for (let j = 1; j <= mesh.Ny; j++) {
      for (let i = 1; i <= mesh.Nx; i++) {
        const idx = j * stride + i;
        mesh.data[idx] =
          this.params.baseValue +
          (this.params.maxVar - this.params.minVar) * Math.random() +
          this.params.minVar;
      }
    }
  }
}

class Mesh {
  constructor(Nx, Ny, dx, dy, data, boundCons) {
    this.Nx = Nx;
    this.Ny = Ny;
    this.dx = dx;
    this.dy = dy;
    this.inv_dxdy = 1 / (dx * dy);
    this.data = data || new Array((Nx + 2) * (Ny + 2)).fill(0);
    this.boundCons = boundCons || new BoundaryConditions();
  }
  load_boundary_conditions(bcs) {
    bcs.forEach((bc) => this.boundCons.addBC(bc));
  }
  init_geometry(initGeom) {
    initGeom.apply(this);
  }
  get(x, y) {
    return this.data[(y + 1) * (this.Nx + 2) + x + 1];
  }
  findFirstNaN() {
    const stride = this.Nx + 2;
    for (let j = 1; j <= this.Ny; j++) {
      for (let i = 1; i <= this.Nx; i++) {
        const idx = j * stride + i;
        if (Number.isNaN(this.data[idx])) {
          return {
            x: i - 1,
            y: j - 1,
            idx,
          };
        }
      }
    }
    return undefined;
  }
  output_mesh(ker_fun) {
    const new_mesh = new Mesh(
      this.Nx,
      this.Ny,
      this.dx,
      this.dy,
      undefined,
      this.boundCons,
    );
    const stride = this.Nx + 2;
    for (let j = 1; j < this.Ny + 1; j++) {
      for (let i = 1; i < this.Nx + 1; i++) {
        const left = j * stride + i - 1;
        const right = j * stride + i + 1;
        const down = (j - 1) * stride + i;
        const up = (j + 1) * stride + i;
        const center = j * stride + i;
        new_mesh.data[center] = ker_fun(
          this.data[left],
          this.data[right],
          this.data[up],
          this.data[down],
          this.data[center],
          this.inv_dxdy,
        );
      }
    }
    new_mesh.boundCons.applyAll(new_mesh);
    return new_mesh;
  }
  write_VTK(istep, folder_name) {
    const vtkContent = `# vtk DataFile Version 3.0
Spinodal Decomposition Step ${istep}
ASCII
DATASET STRUCTURED_GRID
DIMENSIONS ${this.Nx} ${this.Ny} 1
POINTS ${this.Nx * this.Ny} float
${Array.from({ length: this.Ny }, (_, j) =>
  Array.from(
    { length: this.Nx },
    (_, i) => `${i * this.dx} ${j * this.dy} 0`,
  ).join("\n"),
).join("\n")}
POINT_DATA ${this.Nx * this.Ny}
SCALARS CON float
LOOKUP_TABLE default
${Array.from({ length: this.Ny }, (_, j) =>
  Array.from(
    { length: this.Nx },
    (_, i) => `${this.data[(j + 1) * (this.Nx + 2) + i + 1]}`,
  ).join("\n"),
).join("\n")}
`;

    fs.writeFileSync(`${folder_name}/step_${istep}.vtk`, vtkContent);
  }
}

class BulkFreeEnergy {
  constructor(type, params) {
    this.type = type;
    this.params = params;
  }
}

class DoubleWellBulk extends BulkFreeEnergy {
  constructor(A) {
    super("double_well", { A });
  }
  dbulk_dc() {
    return (center_val) =>
      2 * this.params.A * center_val * (1 - center_val) * (1 - 2 * center_val);
  }
}

class InterfacialEnergy {
  constructor(type, params) {
    this.type = type;
    this.params = params;
  }
}

class GradientEnergy extends InterfacialEnergy {
  constructor(kappa) {
    super("gradient", { kappa });
  }
  dint_dc() {
    return (left_val, right_val, up_val, down_val, center_val, inv_dxdy) =>
      -this.params.kappa *
      inv_dxdy *
      (left_val + right_val + up_val + down_val - 4 * center_val);
  }
}

class CH_Solver {
  constructor(mesh, M, dbulk_dc, dint_dc, others = undefined, dt) {
    this.mesh = mesh;
    this.M = M;
    this.dbulk_dc = dbulk_dc;
    this.dint_dc = dint_dc;
    if (others !== undefined) {
      this.others = others;
    }
    this.dt = dt;
  }
  dF_dc() {
    return this.mesh.output_mesh(
      (left_val, right_val, up_val, down_val, center_val, inv_dxdy) => {
        const dbulk = this.dbulk_dc(center_val);
        const dint = this.dint_dc(
          left_val,
          right_val,
          up_val,
          down_val,
          center_val,
          inv_dxdy,
        );
        let dF = dbulk + dint;
        if (this.others !== undefined) {
          dF += this.others(
            left_val,
            right_val,
            up_val,
            down_val,
            center_val,
            inv_dxdy,
          );
        }
        return dF;
      },
    );
  }
  step() {
    const dF = this.dF_dc();
    const dc = dF.output_mesh(
      (left_val, right_val, up_val, down_val, center_val, inv_dxdy) =>
        (left_val + right_val + up_val + down_val - 4 * center_val) * inv_dxdy,
    );
    for (let j = 1; j <= this.mesh.Ny; j++) {
      for (let i = 1; i <= this.mesh.Nx; i++) {
        const idx = j * (this.mesh.Nx + 2) + i;
        this.mesh.data[idx] += this.M * dc.data[idx] * this.dt;
      }
    }
    this.mesh.boundCons.applyAll(this.mesh);
  }
}

export {
  BoundaryConditions,
  CH_Solver,
  DoubleWellBulk,
  GradientEnergy,
  Mesh,
  RandomInit,
};
