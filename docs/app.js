const ramRange = document.querySelector("#ram-range");
const ramNumber = document.querySelector("#ram-number");

const output = {
  memory: document.querySelector("#out-memory"),
  processors: document.querySelector("#out-processors"),
  swap: document.querySelector("#out-swap"),
  default3km: document.querySelector("#out-default"),
  aggressive3km: document.querySelector("#out-aggressive"),
  tone: document.querySelector("#out-tone"),
  config: document.querySelector("#wsl-config-block"),
};

function clamp(value, min, max) {
  return Math.max(min, Math.min(max, value));
}

function recommendationForRam(ram) {
  let windowsHeadroom = 6;
  if (ram >= 64) windowsHeadroom = 12;
  else if (ram >= 32) windowsHeadroom = 8;

  const memory = clamp(ram - windowsHeadroom, 8, Math.max(8, ram - 4));
  const processors = ram >= 64 ? 16 : ram >= 32 ? 12 : ram >= 24 ? 8 : 6;
  const swap = ram <= 16 ? 8 : 16;

  let default3km = "180 to 220";
  let aggressive3km = "not recommended";
  let tone = "Below the pleasant range for big real-data severe-weather WRF. Keep the domain small or use external guidance.";

  if (ram >= 24 && ram < 48) {
    default3km = "220 to 300";
    aggressive3km = "around 300";
    tone = "Usable starter tier. Stay CPU-first and keep output intervals conservative.";
  } else if (ram >= 48 && ram < 96) {
    default3km = "300 to 350";
    aggressive3km = "around 350";
    tone = "This is the sweet spot for a hobbyist 3 km single domain.";
  } else if (ram >= 96) {
    default3km = "350 to 400";
    aggressive3km = "around 400";
    tone = "Aggressive hobbyist tier. Disk and output volume become as important as RAM.";
  }

  return { memory, processors, swap, default3km, aggressive3km, tone };
}

function render() {
  const raw = parseInt(ramNumber.value || ramRange.value || "32", 10);
  const ram = clamp(Number.isFinite(raw) ? raw : 32, 8, 192);
  ramRange.value = ram;
  ramNumber.value = ram;

  const rec = recommendationForRam(ram);

  output.memory.textContent = `${rec.memory} GB`;
  output.processors.textContent = `${rec.processors} threads`;
  output.swap.textContent = `${rec.swap} GB`;
  output.default3km.textContent = rec.default3km;
  output.aggressive3km.textContent = rec.aggressive3km;
  output.tone.textContent = rec.tone;
  output.config.textContent = `[wsl2]
memory=${rec.memory}GB
processors=${rec.processors}
swap=${rec.swap}GB`;
}

function syncFromRange() {
  ramNumber.value = ramRange.value;
  render();
}

function syncFromNumber() {
  ramRange.value = ramNumber.value;
  render();
}

ramRange?.addEventListener("input", syncFromRange);
ramNumber?.addEventListener("input", syncFromNumber);

document.querySelectorAll("[data-jump]").forEach((button) => {
  button.addEventListener("click", () => {
    const target = document.querySelector(button.getAttribute("data-jump"));
    target?.scrollIntoView({ behavior: "smooth", block: "start" });
  });
});

render();
