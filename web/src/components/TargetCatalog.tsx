import { useState } from "react";
import { useQuery } from "@tanstack/react-query";
import { api } from "../api";
import { ErrorMessage, Loading, formatCount } from "./Status";

export function TargetCatalog({
  selectedId,
  onSelect,
}: {
  selectedId: string | null;
  onSelect: (id: string) => void;
}) {
  const [input, setInput] = useState("");
  const [filter, setFilter] = useState({
    query: "",
    changed: false,
    offset: 0,
  });
  const { data, error, isFetching } = useQuery({
    queryKey: ["targets", filter],
    queryFn: ({ signal }) =>
      api.targets(filter.query, filter.changed, filter.offset, signal),
  });
  return (
    <section className="panel catalog">
      <div className="panel-title">
        <span className="step">01</span>
        <h2>Choose a protein</h2>
      </div>
      <form
        onSubmit={(event) => {
          event.preventDefault();
          setFilter({ ...filter, query: input.trim(), offset: 0 });
        }}
      >
        <label htmlFor="query">Gene or transcript prefix</label>
        <div className="search">
          <input
            id="query"
            value={input}
            onChange={(e) => setInput(e.target.value)}
            placeholder="e.g. A1BG"
            autoComplete="off"
          />
          <button>Search</button>
        </div>
        <label className="check">
          <input
            type="checkbox"
            checked={filter.changed}
            onChange={(e) =>
              setFilter({ ...filter, changed: e.target.checked, offset: 0 })
            }
          />
          Non-reference proteins only
        </label>
      </form>
      <ErrorMessage error={error} />
      {isFetching && <Loading>Searching catalog…</Loading>}
      {data && (
        <>
          <p className="muted" role="status">
            {data.total
              ? `${formatCount(data.offset + 1)}–${formatCount(Math.min(data.offset + data.limit, data.total))} of ${formatCount(data.total)} matching records`
              : "No matching records."}
          </p>
          <div className="target-list" aria-label="Matching proteins">
            {data.records.map((row) => (
              <button
                key={row.target_id}
                className={
                  "target" + (selectedId === row.target_id ? " selected" : "")
                }
                aria-pressed={selectedId === row.target_id}
                onClick={() => onSelect(row.target_id)}
              >
                <strong>{row.gene || "Unnamed gene"}</strong>
                <small>{row.transcript}</small>
                <small className="allele">
                  {row.donor} · {formatCount(row.length)} aa ·{" "}
                  {row.changed ? "non-reference" : "reference-matching"}
                </small>
                <small>Record {row.target_id.slice(-8)}</small>
              </button>
            ))}
          </div>
          <div className="pagination">
            <button
              className="secondary"
              disabled={isFetching || !filter.offset}
              onClick={() =>
                setFilter({
                  ...filter,
                  offset: Math.max(0, filter.offset - 25),
                })
              }
            >
              Previous
            </button>
            <button
              className="secondary"
              disabled={isFetching || filter.offset + 25 >= data.total}
              onClick={() =>
                setFilter({ ...filter, offset: filter.offset + 25 })
              }
            >
              Next
            </button>
          </div>
        </>
      )}
    </section>
  );
}
